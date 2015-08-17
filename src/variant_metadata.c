#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include <htslib/vcf.h>

#include "variant_metadata.h"

//{{{ static int get_rowid_callback(void *row_id,
static int get_rowid_callback(void *row_id,
                             int argc,
                             char **argv,
                             char **col_name)
{
    uint32_t *_row_id = (uint32_t *)row_id;

    if ((argc == 1) && (strcmp("rowid",col_name[0]) == 0))
        *_row_id = atoi(argv[0]);

    return 0;
}
//}}}

//{{{ int count_callback(void *count,
int count_callback(void *count,
                          int argc,
                          char **argv,
                          char **col_name)
{
    uint32_t *_count = (uint32_t *)count;
    *_count = *_count + 1;
    return 0;
}
//}}}

//{{{int get_bin_info_callback(void *val,
int get_bin_info_callback(void *info_p,
                          int argc,
                          char **argv,
                          char **col_name)
{
    if ( (argc != 5) ||
         (strcmp(col_name[0], "bin") != 0) ||
         (strcmp(col_name[1], "lo") != 0) ||
         (strcmp(col_name[2], "hi") != 0) ||
         (strcmp(col_name[3], "lo_bin") != 0) ||
         (strcmp(col_name[4], "hi_bin") != 0) ) {
        fprintf(stderr, "FAIL: get_bin_info_callback\n");
        return 1;
    }

    struct bin_info *info = (struct bin_info *)info_p;

    info->bin = atoi(argv[0]);
    info->lo = atof(argv[1]);
    info->hi = atof(argv[2]);
    info->is_lo = atoi(argv[3]);
    info->is_hi = atoi(argv[4]);
    info->was_set = 1;

    return 0;
}
//}}}

//{{{int get_single_int_callback(void *val,
int get_single_int_callback(void *val,
                                   int argc,
                                   char **argv,
                                   char **col_name)
{
    uint32_t *_val = (uint32_t *)val;

    if (argc == 1)
        *_val = atoi(argv[0]);

    return 0;
}
//}}}

//{{{ static int variant_md_info_callback(void *info,
static int variant_md_info_callback(void *info_p,
                                    int argc,
                                    char **argv,
                                    char **col_name)
{
    struct variant_md_info **info = (struct variant_md_info **)info_p;

    *info = (struct variant_md_info *)malloc(sizeof(struct variant_md_info));

    if (argc != 5) {
        fprintf(stderr, "FAIL: variant_md_info_callback\n");
        return 1;
    }

    (*info)->rowid = atoi(argv[0]);
    (*info)->offset = atoll(argv[1]);
    (*info)->bins = atoi(argv[2]);
    (*info)->type = atoi(argv[3]);
    (*info)->source_file = (char *)malloc((strlen(argv[4]) + 1)*sizeof(char));
    strcpy((*info)->source_file, argv[4]);

    return 0;
}
//}}}

//{{{int register_variant_metadata_index(char *bcf_file_name,
int register_variant_metadata_index(char *bcf_file_name,
                                    char *db_file_name,
                                    char *field_name,
                                    uint64_t start_offset,
                                    int type,
                                    uint32_t num_bins)
{

    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(db_file_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        sqlite3_close(db);
        return -1;
    }

    // Check to see if the variant_md_fields table is already there
    char *q_check_table = "SELECT name FROM sqlite_master "
                          "WHERE type='table' AND name='variant_md_fields'";

    uint32_t num_rows = 0;
    rc = sqlite3_exec(db, q_check_table, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return -1;
    }

    if (num_rows == 0) {
        // Create the table
        char *q_create_table = "CREATE TABLE variant_md_fields ("
                               "name TEXT, "
                               "offset INT, "
                               "bins INT, "
                               "type INT, "
                               "source_file TEXT);";

        rc = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            sqlite3_close(db);
            return -1;
        }
    }

    //See if field has been indexed yet
    char *q_check_field = NULL;
    int r = asprintf(&q_check_field,
                     "SELECT * from variant_md_fields WHERE name = '%s'",
                     field_name);
    num_rows = 0;
    rc = sqlite3_exec(db, q_check_field, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return -1;
    }

    free(q_check_field);

    // If so fail
    if (num_rows != 0) {
        fprintf(stderr,
                "The field \"%s\" has already been indexed.\n",
                field_name);

        sqlite3_close(db);
        return -1;
    }

    // If not put in the new row
    char *q_insert_field;
    r = asprintf(&q_insert_field,
                 "INSERT INTO variant_md_fields "
                 "VALUES ('%s',%llu,%d,%d,'%s');",
                 field_name,
                 start_offset,
                 num_bins,
                 type,
                 bcf_file_name);

    rc = sqlite3_exec(db, q_insert_field, NULL, 0, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return -1;
    }

    free(q_insert_field);

    // Return the rowid
    char *q_get_rowid;
    r = asprintf(&q_get_rowid,
                 "SELECT rowid FROM variant_md_fields WHERE name='%s';",
                 field_name);

    uint32_t rowid = 0;
    rc = sqlite3_exec(db, q_get_rowid, get_rowid_callback, &rowid, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return -1;
    }

    sqlite3_close(db);
    return rowid;
}
//}}}

//{{{ uint32_t add_variant_metadata_float_bins(char *db_file_name,
uint32_t add_variant_metadata_float_bins(char *db_file_name,
                                         int rowid,
                                         float *float_bin_range_lo,
                                         float *float_bin_range_hi,
                                         int actual_num_bins,
                                         int less_than_bin,
                                         int greater_than_bin)
{
    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(db_file_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return -1;
    }

    // Check to see if the variant_md_bins table is already there
    char *q_check_table = "SELECT name FROM sqlite_master "
                          "WHERE type='table' AND name='variant_md_bins'";

    uint32_t num_rows = 0;
    rc = sqlite3_exec(db, q_check_table, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return -1;
    }

    if (num_rows == 0) {
        // Create the table
        char *q_create_table = "CREATE TABLE variant_md_bins ( "
                                "md_rowid INT, "
                                "bin INT, "
                                "lo REAL, "
                                "hi REAL, "
                                "lo_bin INT, "
                                "hi_bin INT);";

        rc = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            sqlite3_close(db);
            return -1;
        }
    }

    char *q_insert;
    int r, num_new_rows = 0;

    // add less_than_bin if needed
    if (less_than_bin == 1) {
        // The less than bin has the "lo_bin" field set to 1, and has a hi
        // value equal to the lo value of the first bin
        r = asprintf(&q_insert,
                     "INSERT INTO variant_md_bins VALUES "
                     "(%d, 0, 0.0, %f, 1, 0);",
                     rowid,
                     float_bin_range_lo[0]);

        rc = sqlite3_exec(db, q_insert, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            sqlite3_close(db);
            return -1;
        }
        num_new_rows += 1;
    }
    // add greater_than_bin if needed
     if (greater_than_bin == 1) {
        // The greater than bin has the "hi_bin" field set to 1, and has a lo
        // value equal to the hi value of the last bin
        r = asprintf(&q_insert,
                     "INSERT INTO variant_md_bins VALUES "
                     "(%d, %d, %f, 0.0, 0, 1);",
                     rowid,
                     actual_num_bins + less_than_bin,
                     float_bin_range_hi[actual_num_bins - 1]);

        rc = sqlite3_exec(db, q_insert, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            sqlite3_close(db);
            return -1;
        }
        num_new_rows += 1;
    }
   
    // add range bins
    uint32_t i;
    for (i = 0; i < actual_num_bins; ++i) {
        r = asprintf(&q_insert,
                     "INSERT INTO variant_md_bins VALUES "
                     "(%d, %d, %f, %f, 0, 0);",
                     rowid,
                     i + less_than_bin,
                     float_bin_range_lo[i],
                     float_bin_range_hi[i]);

        rc = sqlite3_exec(db, q_insert, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            sqlite3_close(db);
            return -1;
        }
        num_new_rows += 1;
    }

    free(q_insert);

    rc = sqlite3_close(db);

    return num_new_rows;
}
//}}}

//{{{ int get_variant_metadata_bin_info(char *variant_db_name,
int get_variant_metadata_bin_info(char *variant_db_name,
                                  char *field_name,
                                  uint64_t *offset,
                                  uint32_t *num_bins,
                                  int *type,
                                  char **source_file)
{
    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(variant_db_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return -1;
    }

    char *info_query;

    int r = asprintf(&info_query,
                     "SELECT rowid, offset, bins, type, source_file "
                     "FROM variant_md_fields "
                     "WHERE name == '%s';", field_name);

    struct variant_md_info *info = NULL;

    rc = sqlite3_exec(db,
                      info_query,
                      variant_md_info_callback,
                      &info,
                      &err_msg);

    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        rc = sqlite3_close(db);
        exit(1);
    }

    if (info == NULL) {
        fprintf(stderr,
                "The field \"%s\" has not been indexed.\n",
                field_name);
        rc = sqlite3_close(db);
        return -1;
    }

    *offset = info->offset;
    *num_bins = info->bins;
    *source_file = info->source_file;
    *type = info->type;

    int rowid = info->rowid;

    free(info);
    free(info_query);

    rc = sqlite3_close(db);

    return rowid;
}
//}}}

//{{{ int get_variant_metadata_query_bins(char *variant_db_name,
/*
 * Use an op flag
 * 1: equal
 * 2: less_than
 * 4: less_than_equal
 * 8: greater_than
 * 16: greater_than_equal
 * 32: not_equal
 *
 * lo and hi denote the boundaries of the range.  When the opperation is < (4)
 * or <= (4) are set, then all of the bins that are less than or less than or
 * equal to are returned, and the variable that denotes the upper boundary is 
 * "hi_bin"
 *
 * if 1 (==) or 32 (!=) are set the lo value is used
 * if 2 (<) or 4 (<=) are set lo value is used
 * if 8 (>) or 16 (>=) are set hi value is used
 *
 * In if lo_bin is set the hi value denotes the largest value in the bin,
 * everything smaller than this number is in the lo bin
 *
 * In if hi_bin is set the lo value denotes the lowest value in the bin,
 * everything larger than this number is the in hi bin
 */
int get_variant_metadata_query_bins(char *variant_db_name,
                                    char *field_name,
                                    int rowid,
                                    uint32_t num_bins,
                                    int type,
                                    int op_flag,
                                    float upper_range,
                                    float lower_range,
                                    uint32_t **bins)
{

    /*
     * All of the comparisons will be based on which bins equal the given
     * query.
     * = is the given set
     * < are all bins to the left
     * <= is the give set plus all bins to the left
     * > are all bins to the right
     * >= is the give set plus all bins to the right
     * != are the bins the the left and right
     */

    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(variant_db_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return -1;
    }
    
    // Get matching bins
    /*
     * ==
     * lo 0   1   2   3 hi    lo 0   1   2   3 hi  lo 0   1   2   3 hi
     *    x---o |            |   x---o                x---o           |
     *        x-|-o          |       x---o                x---o       |
     *          | x---o      |           x---o                x---o   |
     *          |            |                                        |   
     *          |            |                                        |  
     *
     */
    char *upper_range_query;
    int r = asprintf(&upper_range_query,
                     "SELECT DISTINCT bin, lo, hi, lo_bin, hi_bin "
                     "FROM variant_md_bins "
                     "WHERE (md_rowid == %d) AND "
                     "( (lo <= %f AND hi > %f) OR "
                     "  (hi > %f AND lo_bin == 1) OR "
                     "  (lo < %f AND hi_bin == 1) )",
                     rowid, upper_range, upper_range, upper_range, upper_range);
    struct bin_info upper_range_info;
    upper_range_info.is_hi = 0;
    upper_range_info.is_lo = 0;
    upper_range_info.was_set = 0;

    rc = sqlite3_exec(db,
                      upper_range_query,
                      get_bin_info_callback,
                      &upper_range_info,
                      &err_msg);

    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }

    char *lower_range_query;
    r = asprintf(&lower_range_query,
                 "SELECT DISTINCT bin, lo, hi, lo_bin, hi_bin "
                 "FROM variant_md_bins "
                 "WHERE (md_rowid == %d) AND "
                 "( (lo <= %f AND hi > %f) OR "
                 "  (hi > %f AND lo_bin == 1) OR "
                 "  (lo < %f AND hi_bin == 1) )",
                 rowid, lower_range, lower_range, lower_range, lower_range);

    struct bin_info lower_range_info;
    lower_range_info.is_hi = 0;
    lower_range_info.is_lo = 0;
    lower_range_info.was_set = 0;

    rc = sqlite3_exec(db,
                      lower_range_query,
                      get_bin_info_callback,
                      &lower_range_info,
                      &err_msg);

    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }

    /*
    * Use an op flag
    * 1: equal
    * 2: less_than
    * 4: less_than_equal
    * 8: greater_than
    * 16: greater_than_equal
    * 32: not_equal
    *
    * if 1 (==) or 32 (!=) are set the lo value is used
    * if 2 (<) or 4 (<=) are set lo value is used
    * if 8 (>) or 16 (>=) are set hi value is used
    */

    // If <, check if we need to step back a bin
    if ( ((op_flag & 2) > 0) &&
         (upper_range_info.is_lo == 0) &&
         (upper_range_info.is_hi == 0) &&
         (upper_range == upper_range_info.lo) ) {
        // Need to move back one bin
        
        upper_range_info.was_set = 0;

        char *step_back_query;
        r = asprintf(&step_back_query,
                     "SELECT DISTINCT bin, lo, hi, lo_bin, hi_bin "
                     "FROM variant_md_bins "
                     "WHERE (md_rowid == %d) AND (bin == %d);",
                     rowid,
                     upper_range_info.bin - 1);;

        rc = sqlite3_exec(db,
                          step_back_query,
                          get_bin_info_callback,
                          &upper_range_info,
                          &err_msg);

        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            exit(1);
        }

        free(step_back_query);
    }

    // If >, check if we need to step forward a bin
    if ( ((op_flag & 8) > 0) &&
         ((type == BCF_HT_INT) || (type == BCF_HT_FLAG)) &&
         (lower_range_info.is_lo == 0) &&
         (lower_range_info.is_hi == 0) &&
         (lower_range == lower_range_info.lo) &&
         (( ((int)lower_range_info.hi) - ((int)lower_range_info.lo)) == 1) ) {
        // Need to move forward one bin

        lower_range_info.was_set = 0;

        char *step_forward_query;
        r = asprintf(&step_forward_query,
                     "SELECT DISTINCT bin, lo, hi, lo_bin, hi_bin "
                     "FROM variant_md_bins "
                     "WHERE (md_rowid == %d) AND (bin == %d);",
                     rowid,
                     lower_range_info.bin + 1);;

        rc = sqlite3_exec(db,
                          step_forward_query,
                          get_bin_info_callback,
                          &lower_range_info,
                          &err_msg);

        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            sqlite3_free(err_msg);
            exit(1);
        }

        free(step_forward_query);
    }



    uint32_t i,bin_i = 0;
    if (op_flag == 1) { // == 

        // Check to see if no bins were found
        if (upper_range_info.was_set == 0) {
            if ( (type == BCF_HT_INT) || (type == BCF_HT_FLAG))
                fprintf(stderr, "Range for %s == %d is EMPTY.\n",
                        field_name,
                        (int)upper_range);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s == %f is EMPTY.\n",
                        field_name,
                        upper_range);
 
            return 0;
        }

        if (type == BCF_HT_INT)
            fprintf(stderr, "Range for %s == %d is [%d, %d).\n",
                    field_name,
                    (int)upper_range,
                    (int)upper_range_info.lo,
                    (int)upper_range_info.hi);
        else if (type == BCF_HT_REAL)
            fprintf(stderr, "Range for %s == %f is [%f, %f).\n",
                    field_name,
                    upper_range,
                    upper_range_info.lo,
                    upper_range_info.hi);

        *bins = (uint32_t *)malloc(sizeof(uint32_t));
        (*bins)[0] = upper_range_info.bin;
        return 1;

    } else if ((op_flag == 2) || (op_flag ==4)) { // < <=
        char *op = "=";
        if (op_flag == 2)
            op = "";

        if (type == BCF_HT_FLAG) {
                fprintf(stderr, 
                        ">,<=,>, and > not defined for FLAG data type.\n");
                return 0;
        }

        // Check to see if no bins were found
        if (upper_range_info.was_set == 0) {
            //if ( (type == BCF_HT_INT) || (type == BCF_HT_FLAG))
            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s <%s %d is EMPTY.\n",
                        field_name,
                        op,
                        (int)upper_range);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s <%s %f is EMPTY.\n",
                        field_name,
                        op,
                        upper_range);
            return 0;
        }

        if (upper_range_info.is_lo == 1) {

            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s <%s %d is (-inf, %d).\n",
                        field_name,
                        op,
                        (int)upper_range,
                        (int)upper_range_info.hi);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s <%s %f is (-inf, %f).\n",
                        field_name,
                        op,
                        upper_range,
                        upper_range_info.hi);

            *bins = (uint32_t *)malloc(sizeof(uint32_t));
            (*bins)[0] = upper_range_info.bin;
            return 1;
        } else if (upper_range_info.is_hi == 1) {

            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s <%s %d is -inf ... (%d,inf).\n",
                        field_name,
                        op,
                        (int)upper_range,
                        (int)upper_range_info.lo);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s <%s %f is -inf ... (%f, inf).\n",
                        field_name,
                        op,
                        upper_range,
                        upper_range_info.lo);


            *bins = (uint32_t *)malloc(num_bins*sizeof(uint32_t));
            for (i = 0; i < num_bins; ++i)
                (*bins)[i] = i;
            return num_bins;
        } else {

            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s <%s %d is -inf ... [%d, %d).\n",
                        field_name,
                        op,
                        (int)upper_range,
                        (int)upper_range_info.lo,
                        (int)upper_range_info.hi);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s <%s %f is -inf ... [%f, %f).\n",
                        field_name,
                        op,
                        upper_range,
                        upper_range_info.lo,
                        upper_range_info.hi);

            *bins = (uint32_t *)malloc(upper_range_info.bin+1*sizeof(uint32_t));
            for (i = 0; i <= upper_range_info.bin; ++i) {
                (*bins)[bin_i] = i;
                bin_i+=1;
            }
            return bin_i;
        }
    } else if ( (op_flag == 8) || (op_flag == 16)) { // >
        char *op = "=";
        if (op_flag == 8)
            op = "";

        if (type == BCF_HT_FLAG) {
                fprintf(stderr, 
                        ">,<=,>, and > not defined for FLAG data type.\n");
                return 0;
        }

        // Check to see if no bins were found
        if (lower_range_info.was_set == 0) {
            //if ( (type == BCF_HT_INT) || (type == BCF_HT_FLAG))
            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s >%s %d is EMPTY.\n",
                        field_name,
                        op,
                        (int)lower_range);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s >%s %f is EMPTY.\n",
                        field_name,
                        op,
                        lower_range);
            return 0;
        }

        if (lower_range_info.is_hi == 1) {

            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s >%s %d is (%d, inf).\n",
                        field_name,
                        op,
                        (int)lower_range,
                        (int)lower_range_info.lo);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s >%s %f is (%f, inf).\n",
                        field_name,
                        op,
                        lower_range,
                        lower_range_info.lo);

            *bins = (uint32_t *)malloc(sizeof(uint32_t));
            (*bins)[0] = lower_range_info.bin;
            return 1;
        } else if (lower_range_info.is_lo == 1) {

            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s >%s %d is (-inf, %d) ... inf.\n",
                        field_name,
                        op,
                        (int)lower_range,
                        (int)lower_range_info.hi);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s >%s %f is (-inf, %f) ... inf.\n",
                        field_name,
                        op,
                        lower_range,
                        lower_range_info.hi);

            *bins = (uint32_t *)malloc(num_bins*sizeof(uint32_t));
            for (i = 0; i < num_bins; ++i)
                (*bins)[i] = i;
            return num_bins;
        } else {

            //if ((type == BCF_HT_INT) || (type == BCF_HT_FLAG))
            if (type == BCF_HT_INT)
                fprintf(stderr, "Range for %s >%s %d is [%d, %d) ... inf.\n",
                        field_name,
                        op,
                        (int)lower_range,
                        (int)lower_range_info.lo,
                        (int)lower_range_info.hi);
            else if (type == BCF_HT_REAL)
                fprintf(stderr, "Range for %s >%s %f is [%f, %f) ... inf.\n",
                        field_name,
                        op,
                        lower_range,
                        lower_range_info.lo,
                        lower_range_info.hi);

            *bins = (uint32_t *)malloc(
                    (num_bins-lower_range_info.bin)*sizeof(uint32_t));
            for (i = lower_range_info.bin; i < num_bins; ++i) {
                (*bins)[bin_i] = i;
                bin_i+=1;
            }
            return bin_i;
        }
    } else if (op_flag == 32) { // != 

        if ( (type == BCF_HT_INT) || (type == BCF_HT_FLAG))
            fprintf(stderr, "Range for %s != %d is NOT [%d, %d).\n",
                    field_name,
                    (int)upper_range,
                    (int)upper_range_info.lo,
                    (int)upper_range_info.hi);
        else if (type == BCF_HT_REAL)
            fprintf(stderr, "Range for %s != %f is NOT [%f, %f).\n",
                    field_name,
                    upper_range,
                    upper_range_info.lo,
                    upper_range_info.hi);

        *bins = (uint32_t *)malloc((num_bins - 1)*sizeof(uint32_t));
        for (i = 0; i < upper_range_info.bin; ++i) {
            (*bins)[bin_i] = i;
            bin_i+=1;
        }
        for (i = upper_range_info.bin+1; i < num_bins; ++i) {
            (*bins)[bin_i] = i;
            bin_i+=1;
        }
        return bin_i;
    } else if ( ((op_flag & 2+4) > 0) && ((op_flag & 8+16) > 0) ) {

        char *upper_op = "=";
        if ((op_flag & 2) > 0)
            upper_op = "";
        char *lower_op = "=";
        if ((op_flag & 8) > 0)
            lower_op = "";


        if (type == BCF_HT_FLAG) {
                fprintf(stderr, 
                        ">,<=,>, and > not defined for FLAG data type.\n");
                return 0;
        }


        //if ((type == BCF_HT_INT) || (type == BCF_HT_FLAG))
        if (type == BCF_HT_INT)
            fprintf(stderr, "Range for %d <%s %s <%s %d is ",
                    (int)lower_range,
                    lower_op,
                    field_name,
                    upper_op,
                    (int)upper_range);
        else if (type == BCF_HT_REAL)
            fprintf(stderr, "Range for %f <%s  %s <%s %f is ",
                    lower_range,
                    lower_op,
                    field_name,
                    upper_op,
                    upper_range);

        if ( (upper_range_info.bin-lower_range_info.bin+1) <= 0) {
            fprintf(stderr, "EMPTY\n");
            return 0;
        }
        
        if (lower_range_info.is_lo == 1) {
            fprintf(stderr, "-inf ...");
        } else if (lower_range_info.is_hi == 1) {
            if (type == BCF_HT_INT) {
                fprintf(stderr,
                        "(%d, inf) ... ",
                        (int)lower_range_info.lo);
            } else if (type == BCF_HT_REAL) {
                fprintf(stderr,
                        "(%f, inf) ... ",
                        lower_range_info.lo);
            }
        } else {
            if (type == BCF_HT_INT) {
                fprintf(stderr,
                        "(%d, %d) ... ",
                        (int)lower_range_info.lo,
                        (int)lower_range_info.hi);
            } else if (type == BCF_HT_REAL) {
                fprintf(stderr,
                        "(%f, %f) ... ",
                        lower_range_info.lo,
                        lower_range_info.hi);
            }
        }

        if (upper_range_info.is_hi == 1) {
            fprintf(stderr,"inf\n");
        }else if (upper_range_info.is_lo == 1) {
            if (type == BCF_HT_INT) {
                fprintf(stderr,
                        "(-inf, %d)\n",
                        (int)upper_range_info.hi);
            } else if (type == BCF_HT_REAL) {
                fprintf(stderr,
                        "(inf, %f)\n",
                        upper_range_info.hi);
            }
        } else {
            if (type == BCF_HT_INT) {
                fprintf(stderr,
                        "(%d, %d)\n",
                        (int)upper_range_info.lo,
                        (int)upper_range_info.hi);
            } else if (type == BCF_HT_REAL) {
                fprintf(stderr,
                        "(%f, %f)\n",
                        upper_range_info.lo,
                        upper_range_info.hi);
            }
        }

        *bins = (uint32_t *)malloc(
                (upper_range_info.bin-lower_range_info.bin+1)*sizeof(uint32_t));
        for (i = lower_range_info.bin; i <= upper_range_info.bin; ++i) {
            (*bins)[bin_i] = i;
            bin_i+=1;
        }
        return bin_i;
    }

    sqlite3_close(db);
    return 0;
}
//}}}
//
//
//TODO THIS
#if 0
//{{{ uint32_t get_wah_bitmap_in_place(struct wah_file wf,
uint32_t get_wah_bitmap_inline_index_in_place(FILE *f, 
                                              uint64_t offset,
                                              uint64_t num_bitmaps,
                                              uint32_t bitmap_id,
                                              uint32_t **wah_sizes,
                                              uint32_t *wah_bitmap)
{

    // get the size of the WAH-encoded bitmap
    uint64_t wah_size = 0;
    uint64_t wah_offset = 0;
    if ((wah_record == 0) && (bitmap == 0)) {
        wah_size = wf.record_offsets[wah_record + bitmap];
        wah_offset = wf.header_offset;
    } else {
        wah_size = wf.record_offsets[wah_record*4 + bitmap] - 
                   wf.record_offsets[wah_record*4 + bitmap - 1];

        wah_offset = wf.header_offset +
                     sizeof(uint32_t) * 
                        (wf.record_offsets[wah_record*4 + bitmap] - wah_size);
    }

    fseek(wf.file, wah_offset, SEEK_SET);
    int r = fread(*wah_bitmap,sizeof(uint32_t),wah_size,wf.file);

    return (uint32_t)wah_size;
}
//}}}
#endif
