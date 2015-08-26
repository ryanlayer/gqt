#define _GNU_SOURCE
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include <htslib/vcf.h>
#include "ped.h"

//{{{ static int uint32_t_ll_callback(void *ll_p,
static int uint32_t_ll_callback(void *ll_p,
                    int argc,
                    char **argv,
                    char **col_name)
{
    struct uint32_t_ll *ll = (struct uint32_t_ll *)ll_p;

    struct uint32_t_ll_node *new_node = (struct uint32_t_ll_node *)
        malloc(sizeof(struct uint32_t_ll_node));
    if (!new_node )
        err(EX_OSERR, "malloc error");
    new_node->v = atoi(argv[0]);
    new_node->next = NULL;

    if (ll->head == NULL)
        ll->head = new_node;
    else
        ll->tail->next = new_node;

    ll->tail = new_node;

    ll->len = ll->len + 1;
    return 0;
}
//}}}

//{{{ static int char_ll_callback(void *ll_p,
static int char_ll_callback(void *ll_p,
                    int argc,
                    char **argv,
                    char **col_name)
{
    if (argc != 2) {
        fprintf(stderr,
                "FAILURE: Cannot get label value. Expecte 2 columns, "
                "but recieved %d.\n",
                argc);
        exit(1);
    }
    struct char_ll *ll = (struct char_ll *)ll_p;

    struct char_ll_node *new_node = (struct char_ll_node *)
        malloc(sizeof(struct char_ll_node));
    if (!new_node )
        err(EX_OSERR, "malloc error");

    uint32_t bcf_i = 0;
    uint32_t label_i = 1;

    if (strlen(argv[label_i]) == 0) {
        fprintf(stderr,
                "FAILURE: Blank label found for BCF_ID:%s\n",
                argv[label_i]);
        exit(1);
    }

    new_node->v = (char *) malloc(strlen(argv[label_i])*sizeof(char));
    if (!new_node->v )
        err(EX_OSERR, "malloc error");
    strcpy(new_node->v, argv[label_i]);
    new_node->next = NULL;

    if (ll->head == NULL)
        ll->head = new_node;
    else
        ll->tail->next = new_node;

    ll->tail = new_node;

    ll->len = ll->len + 1;
    return 0;
}
//}}}

//{{{ uint32_t convert_file_by_name_ped_to_db(char *bcf_file_name,
/*
 * Every sample DB will have BCF_Sample and BCF_ID fields that are defined by
 * the BCF. Extra fields can be included by the ped_file_name.
 *
 */
uint32_t convert_file_by_name_ped_to_db(char *bcf_file_name,
                                        char *ped_file_name,
                                        uint32_t col,
                                        char *db_name)
{

    char **ped_field_names;
    int *ped_field_types;
    uint32_t i, j, num_ped_fields = 0;

    // Figure out which fields will be in the DB
    if (ped_file_name != NULL) {

        FILE *ped_f = fopen(ped_file_name, "r");
        if (!ped_f)
            err(EX_NOINPUT, "Cannot open file '%s'", ped_file_name);

        char *line = NULL, *tmp_line;
        size_t len = 0;

        ssize_t read = getline(&line, &len, ped_f);
        if (read == -1) {
            if (feof(ped_f))
                errx(EX_NOINPUT,
                     "Error reading file '%s': End of file",
                     ped_file_name);
            err(EX_NOINPUT, "Error reading file '%s'", ped_file_name);
        }

        // Scan the first line to get the names of the fields
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        if (line[0] == '#'){
            line++;
        }

        char *word;
        tmp_line = (char *) malloc(strlen(line) * sizeof(char));
        if (!tmp_line )
            err(EX_OSERR, "malloc error");
        strcpy(tmp_line, line);

        word = strtok(tmp_line, "\t");
        while (word != NULL) {
            num_ped_fields += 1;
            word = strtok(NULL, "\t");
        }

        if (num_ped_fields == 0) {
            errx(EX_NOINPUT, "Empty PED file '%s'.", ped_file_name);
        } else if (num_ped_fields < col) {
            errx(EX_NOINPUT,
                 "Too few columns in PED file '%s'. Sample IDs column "
                 "specified as %d,  but only %d columns present.",
                 ped_file_name,
                 col,
                 num_ped_fields);
        } else {
            ped_field_names = (char **) malloc(num_ped_fields * sizeof(char *));
            if (!ped_field_names )
                err(EX_OSERR, "malloc error");
            strcpy(tmp_line, line);

            // Set field names
            ped_field_names[0] = strtok(tmp_line, "\t");
            for (i = 1; i < num_ped_fields; ++i)
                ped_field_names[i] = strtok(NULL, "\t");

            // Convert " " to "_"
            for (i = 0; i < num_ped_fields; ++i) {
                for (j = 0; j < strlen(ped_field_names[i]); ++j) {
                    if (ped_field_names[i][j] == ' ')
                        ped_field_names[i][j] = '_';
                }
            }

            // Check for problems with field names
            for (i = 0; i < num_ped_fields; ++i) {
                int r = check_field_name(ped_field_names[i]);
                if (r >= 0) {
                    errx(EX_NOINPUT, 
                         "Invalid character '%c' in field name '%s' from file "
                         "'%s'",
                         ped_field_names[i][r],
                         ped_field_names[i],
                         ped_file_name);
                }
            }

            // Set field types
            ped_field_types = (int *) malloc(num_ped_fields * sizeof(int));
            if (!ped_field_types )
                err(EX_OSERR, "malloc error");

            for (i = 0; i < num_ped_fields; ++i)
                ped_field_types[i] = 1;

            uint32_t line_no = 2;

            while ( (read = getline(&line, &len, ped_f)) != -1) {
                if (line[strlen(line) - 1] == '\n')
                    line[strlen(line) - 1] = '\0';

                uint32_t j;
                word = strtok(line, "\t");

                for (i = 0; i < num_ped_fields; ++i) {
                    if (word == NULL) {
                        errx(EX_NOINPUT,
                             "Missing field in file '%s' on line %u.\n",
                             ped_file_name,
                             line_no);
                    }

                    int v;
                    ped_field_types[i] &= is_int(word, &v);
                    word = strtok(NULL, "\t");
                }

                // unparsed data on this line
                if (word != NULL)
                    errx(EX_NOINPUT,
                         "Extra field in file '%s' on line %u.\n",
                         ped_file_name,
                         line_no);

                line_no += 1;
            }

            // check to see if no data is read
            if (line_no == 2) 
                errx(EX_NOINPUT, "No data in PED file '%s'", ped_file_name);

            if (ferror(ped_f) != 0 )
                err(EX_NOINPUT, "Error reading file '%s'", ped_file_name);
        }
        fclose(ped_f);
    }

    fprintf(stderr, "Creating sample database %s\n", ped_file_name);

    if (num_ped_fields > 0 ) {
        fprintf(stderr,
                "Adding the following fields from %s\n",
                ped_file_name);


        for (i = 0; i < num_ped_fields; ++i)
            fprintf(stderr,
                    "%s\t%s\n",
                    ped_field_names[i],
                    ped_field_types[i] ? "INT" : "TEXT");
    }


    fprintf(stderr, "Adding the following fields from %s\n",
                bcf_file_name);
    fprintf(stderr, "BCF_ID\tINT\nBCF_Sample\tTEXT\n");


    // Add the minimum fields
    char *q_create_table, *q_create_table_tmp;
    int r = asprintf(&q_create_table,
                     "CREATE TABLE ped(BCF_ID INTEGER, BCF_Sample TEXT");
    if (r == -1) err(EX_OSERR, "asprintf error");

    // Add fields from PED
    for (i = 0; i < num_ped_fields; ++i) {
        if (ped_field_types[i] == 1) {
            r = asprintf(&q_create_table_tmp,
                         "%s, %s INTEGER",
                         q_create_table,
                         ped_field_names[i]);

            if (r == -1) err(EX_OSERR, "asprintf error");

        } else {
            r = asprintf(&q_create_table_tmp,
                         "%s, %s TEXT",
                         q_create_table,
                         ped_field_names[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");
        }

        free(q_create_table);
        q_create_table = q_create_table_tmp;
    }

    r = asprintf(&q_create_table_tmp, "%s);", q_create_table);
    if (r == -1) err(EX_OSERR, "asprintf error");
    free(q_create_table);
    q_create_table = q_create_table_tmp;

    // removed DB if it is there
    struct stat buffer;
    r = stat(db_name, &buffer);

    if (r == 0)
        remove(db_name);

    // create the table
    sqlite3 *db;
    char *err_msg = NULL;
    int rc = sqlite3_open(db_name, &db);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,
            "SQL error '%s' for database '%s'",
            err_msg, db_name);


    rc = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,
            "SQL error '%s' in query '%s'",
            err_msg,
            q_create_table);

    // open BCF
    htsFile *fp = hts_open(bcf_file_name,"rb");
    if ( !fp ) 
        err(EX_DATAERR, "Could not read file: %s", bcf_file_name);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if ( !hdr )
        err(EX_DATAERR, "Could not read the header: %s", bcf_file_name);

    // Add the sample names and location from the BCF file
    char *q;

    char *base_insert_txt = "INSERT INTO ped(BCF_ID, BCF_Sample)"
                            "VALUES (?, ?);";
    sqlite3_stmt *base_insert_stmt = NULL;
    for (i = 0; i < hdr->n[BCF_DT_SAMPLE]; ++i) {
        rc = sqlite3_prepare_v2(db,
                                base_insert_txt,
                                strlen(base_insert_txt),
                                &base_insert_stmt,
                                NULL);
        if (rc != SQLITE_OK) 
            errx(EX_SOFTWARE,
                "Can't prepare insert statment %s (%i): %s",
                base_insert_txt, rc, sqlite3_errmsg(db));

        rc = sqlite3_bind_int(base_insert_stmt, 1, i);
        if (rc != SQLITE_OK) 
            errx(EX_SOFTWARE,
                "Error binding value in insert (%i): %s",
                rc, sqlite3_errmsg(db));


        rc = sqlite3_bind_text(base_insert_stmt,
                               2,
                               hdr->samples[i],
                               strlen(hdr->samples[i]),
                               NULL);
        if (rc != SQLITE_OK) 
            errx(EX_SOFTWARE,
                "Error binding value in insert (%i): %s",
                rc, sqlite3_errmsg(db));

        rc = sqlite3_step(base_insert_stmt);
        if (rc != SQLITE_DONE)
            errx(EX_SOFTWARE,
                "Error on insert(%i): %s",
                rc, sqlite3_errmsg(db));

        rc = sqlite3_finalize(base_insert_stmt);
        if (rc != SQLITE_OK)
            errx(EX_SOFTWARE,
                "Error finalizing insert(%i): %s",
                rc, sqlite3_errmsg(db));
    }

    if (num_ped_fields > 0) {
        fprintf(stderr,
                "Joining values based on BCF_Sample in %s and %s in %s.\n",
                bcf_file_name,
                ped_field_names[col - 1],
                ped_file_name);

        FILE *ped_f = fopen(ped_file_name, "r");
        if (!ped_f)
            err(EX_NOINPUT, "Cannot open file '%s'", ped_file_name);

        char *line = NULL;
        size_t len = 0;

        ssize_t read = getline(&line, &len, ped_f); // skip header
        if (read == -1) {
            if (feof(ped_f))
                errx(EX_NOINPUT,
                     "Error reading file '%s': End of file",
                     ped_file_name);
            err(EX_NOINPUT, "Error reading file '%s'", ped_file_name);
        }

        char *update_query_txt = "UPDATE ped SET ";
        for (i = 0; i < num_ped_fields; ++i) {
            if (i != 0) {
                r = asprintf(&update_query_txt,"%s,",update_query_txt);
                if (r == -1) err(EX_OSERR, "asprintf error");
            }

            r = asprintf(&update_query_txt,
                         "%s %s=?",
                         update_query_txt,
                         ped_field_names[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");
        }

        r = asprintf(&update_query_txt,
                     "%s WHERE BCF_Sample == ?;",
                     update_query_txt);

        if (r == -1) err(EX_OSERR, "asprintf error");

        sqlite3_stmt *update_stmt = NULL;

        char **ped_values = (char **) malloc(num_ped_fields *sizeof(char *));
        if (!ped_values)
            err(EX_OSERR, "malloc error");

        int total_changes = 0, total_lines = 0;

        while ( (read = getline(&line, &len, ped_f)) != -1) {
            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';

            rc = sqlite3_prepare_v2(db,
                                    update_query_txt,
                                    strlen(update_query_txt),
                                    &update_stmt,
                                    NULL);

            if (rc != SQLITE_OK) 
                errx(EX_SOFTWARE,
                    "Can't prepare insert statment %s (%i): %s",
                    update_query_txt, rc, sqlite3_errmsg(db));

            ped_values[0] = strtok(line, "\t");
            for (i = 1; i < num_ped_fields; ++i)
                ped_values[i] = strtok(NULL, "\t");

            for (i = 0; i < num_ped_fields; ++i) {
                if (ped_field_types[i] == 1) {
                    int v;
                    rc = is_int(ped_values[i], &v);

                    // This should never happen
                    if (rc == 0)
                        errx(EX_SOFTWARE,
                             "Value '%s' was erroneously identified as an "
                             "int in %s.",
                             ped_values[i],
                             ped_file_name);

                    rc = sqlite3_bind_int(update_stmt, i+1, v);
                    if (rc != SQLITE_OK) 
                        errx(EX_SOFTWARE,
                             "Error binding value in insert (%i): %s",
                             rc, sqlite3_errmsg(db));

                } else {
                    rc = sqlite3_bind_text(update_stmt,
                                           i+1,
                                           ped_values[i],
                                           strlen(ped_values[i]),
                                           NULL);
                    if (rc != SQLITE_OK) 
                        errx(EX_SOFTWARE,
                             "Error binding value in insert (%i): %s",
                             rc, sqlite3_errmsg(db));
                }
            }

            rc = sqlite3_bind_text(update_stmt,
                                   num_ped_fields + 1,
                                   ped_values[col - 1],
                                   strlen(ped_values[col - 1]),
                                   NULL);
            if (rc != SQLITE_OK) 
                errx(EX_SOFTWARE,
                     "Error binding value in insert (%i): %s",
                     rc, sqlite3_errmsg(db));

            rc = sqlite3_step(update_stmt);
            if (rc != SQLITE_DONE)
                errx(EX_SOFTWARE,
                     "Error on insert(%i): %s",
                     rc, sqlite3_errmsg(db));

            int changes = sqlite3_changes(db);

            total_changes += changes;

            if (changes > 1)
                errx(EX_SOFTWARE,
                     "Multiple rows (%d) changed by update query from "
                     "'%s' and '%s'.",
                     changes,
                     bcf_file_name,
                     ped_file_name);

            if (changes == 0)
                warnx("WARNING: No match found for sample '%s' from PED file.",
                      ped_values[col - 1]);

            rc = sqlite3_finalize(update_stmt);
            if (rc != SQLITE_OK)
                errx(EX_SOFTWARE,
                     "Error finalizing insert(%i): %s",
                     rc, sqlite3_errmsg(db));

            total_lines += 1;
        }

        if (total_changes == 0)
            warnx("WARNING: None of the samples names from column %d in PED "
                  "file %s matched sample names in VCF/BCF %s'",
                  col,
                  ped_file_name,
                  bcf_file_name);

        fprintf(stderr,
                "%d of %d PED samples matched VCF/BCF database records.\n",
                total_changes, total_lines);
        fclose(ped_f);
    }

    sqlite3_close(db);

    return 0;
}
//}}}

//{{{ uint32_t resolve_ind_query(uint32_t **R, char *query, char *ped_db_file)
uint32_t resolve_ind_query(uint32_t **R, char *query, char *ped_db_file)
{
    sqlite3 *db;
    char *err_msg = NULL;
    int rc = sqlite3_open(ped_db_file, &db);
    if( rc != SQLITE_OK )
        err(EX_NOINPUT,
            "SQL error '%s' for database '%s'",
            err_msg, ped_db_file);

    char *test_q;
    int r;
    if (strlen(query) == 0) {
        r = asprintf(&test_q, "SELECT BCF_ID FROM ped ORDER BY BCF_ID");
        if (r == -1) err(EX_OSERR, "asprintf error");
    } else {
        r = asprintf(&test_q,
                     "SELECT BCF_ID FROM ped WHERE %s ORDER BY BCF_ID;",
                     query);
        if (r == -1) err(EX_OSERR, "asprintf error");
    }

    //PRAGMA database.table_info(table-name);

    struct uint32_t_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    rc = sqlite3_exec(db, test_q, uint32_t_ll_callback, &ll, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,"SQL error '%s' in query '%s'", err_msg, test_q);

    *R = (uint32_t *) malloc(ll.len * sizeof(uint32_t));
    if (!R)
        err(EX_OSERR, "malloc error");

    struct uint32_t_ll_node *tmp, *curr = ll.head;
    uint32_t i;
    for (i = 0; i < ll.len; ++i) {
        (*R)[i] = curr->v;
        tmp = curr->next;
        free(curr);
        curr = tmp;
    }

    sqlite3_close(db);
    return ll.len;
}
//}}}

//{{{ uint32_t resolve_label_query(char ***R,
uint32_t resolve_label_query(char ***R,
                             char *label_id,
                             char *query,
                             char *ped_db_file)
{
    sqlite3 *db;
    char *err_msg = NULL;
    int rc = sqlite3_open(ped_db_file, &db);
    if( rc != SQLITE_OK )
        err(EX_NOINPUT,
            "SQL error '%s' for database '%s'",
            err_msg, ped_db_file);

    char *test_q;
    int r;
    if (strlen(query) == 0) {
        r = asprintf(&test_q, "SELECT BCF_ID,%s FROM ped ORDER BY BCF_ID",
                label_id);
        if (r == -1) err(EX_OSERR, "asprintf error");
    } else {
        r = asprintf(&test_q,
                     "SELECT BCF_ID,%s FROM ped WHERE %s ORDER BY BCF_ID;",
                     label_id,
                     query);
        if (r == -1) err(EX_OSERR, "asprintf error");
    }

    struct char_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    rc = sqlite3_exec(db, test_q, char_ll_callback, &ll, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,"SQL error '%s' in query '%s'", err_msg, test_q);

    *R = (char **) malloc(ll.len * sizeof(char *));
    if (!R)
        err(EX_OSERR, "malloc error");

    struct char_ll_node *tmp, *curr = ll.head;
    uint32_t i;
    for (i = 0; i < ll.len; ++i) {
        (*R)[i] = curr->v;
        tmp = curr->next;
        free(curr);
        curr = tmp;
    }

    sqlite3_close(db);
    return ll.len;
}
//}}}
