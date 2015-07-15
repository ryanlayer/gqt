#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include "genotq.h"

//{{{ static int uint32_t_ll_callback(void *ll_p,
static int uint32_t_ll_callback(void *ll_p,
                    int argc,
                    char **argv,
                    char **col_name)
{
    struct uint32_t_ll *ll = (struct uint32_t_ll *)ll_p;

    struct uint32_t_ll_node *new_node = (struct uint32_t_ll_node *)
        malloc(sizeof(struct uint32_t_ll_node));
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

    uint32_t bcf_i = 0;
    uint32_t label_i = 1;

    if (strlen(argv[label_i]) == 0) {
        fprintf(stderr,
                "FAILURE: Blank label found for BCF_ID:%s\n",
                argv[label_i]);
        exit(1);
    }

    new_node->v = (char *) malloc(strlen(argv[label_i])*sizeof(char));
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
            err(EX_NOINPUT, "Cannot open file \"%s\"", ped_file_name);

        fprintf(stderr, "Creating sample database %s\n", ped_file_name);

        fprintf(stderr, "Adding the following fields from %s\n",
                ped_file_name);

        char *line = NULL, *tmp_line;
        size_t len = 0;

        ssize_t read = getline(&line, &len, ped_f);
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        if (line[0] == '#'){
            line++;
        }
        char *word;
        tmp_line = (char *) malloc(strlen(line) * sizeof(char));
        strcpy(tmp_line, line);

        word = strtok(tmp_line, "\t");
        while (word != NULL) {
            num_ped_fields += 1;
            word = strtok(NULL, "\t");
        }

        if (num_ped_fields > 0) {
            ped_field_names = (char **) malloc(num_ped_fields * sizeof(char *));
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

            // Set field types
            ped_field_types = (int *) malloc(num_ped_fields * sizeof(int));
            for (i = 0; i < num_ped_fields; ++i)
                ped_field_types[i] = 1;

            uint32_t line_no = 2;
            while ( (read = getline(&line, &len, ped_f)) != -1) {
                if (line[strlen(line) - 1] == '\n')
                    line[strlen(line) - 1] = '\0';

                uint32_t j;
                //fprintf(stderr, "%s\n", line);
                word = strtok(line, "\t");
                for (i = 0; i < num_ped_fields; ++i) {
                    //fprintf(stderr, "%s\n", word);
                    if (word == NULL) {
                        fprintf(stderr,
                                "ERROR: Missing field in file %s on line %u\n",
                                ped_file_name,
                                line_no);
                        exit(1);
                    }
                    for (j = 0; j < strlen(word); ++j)
                        ped_field_types[i] &= isdigit((int)word[j]);
                    word = strtok(NULL, "\t");
                }
                line_no += 1;
            }
        }
        fclose(ped_f);
    }

    for (i = 0; i < num_ped_fields; ++i)
        fprintf(stderr,
                "%s\t%s\n",
                ped_field_names[i],
                ped_field_types[i] ? "INT" : "TEXT");


    fprintf(stderr, "Adding the following fields from %s\n",
                bcf_file_name);
    fprintf(stderr, "BCF_ID\tINT\nBCF_Sample\tTEXT\n");


    // Add the minimum fields
    char *q_create_table, *q_create_table_tmp;
    int r = asprintf(&q_create_table,
                     "CREATE TABLE ped(BCF_ID INTEGER, BCF_Sample TEXT");

    // Add fields from PED
    for (i = 0; i < num_ped_fields; ++i) {
        if (ped_field_types[i] == 1)
            r = asprintf(&q_create_table_tmp,
                         "%s, %s INTEGER",
                         q_create_table,
                         ped_field_names[i]);

        else
            r = asprintf(&q_create_table_tmp,
                         "%s, %s TEXT",
                         q_create_table,
                         ped_field_names[i]);

        free(q_create_table);
        q_create_table = q_create_table_tmp;
    }

    r = asprintf(&q_create_table_tmp, "%s);", q_create_table);
    free(q_create_table);
    q_create_table = q_create_table_tmp;

    // removed DB if it is there
    struct stat buffer;
    r = stat(db_name, &buffer);

    if (r == 0)
        remove(db_name);

    // create the table
    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(db_name, &db);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,
            "SQL error \"%s\" for database \"%s\"",
            err_msg, db_name);


    rc = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,
            "SQL error \"%s\" in query \"%s\"",
            err_msg,
            q_create_table);


    // open BCF
    htsFile *fp    = hts_open(bcf_file_name,"rb");
    if ( !fp ) 
        err(EX_DATAERR, "Could not read file: %s", bcf_file_name);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if ( !hdr )
        err(EX_DATAERR, "Could not read the header: %s", bcf_file_name);

    // Add the sample names and location from the BCF file
    char *q;
    for (i = 0; i < hdr->n[BCF_DT_SAMPLE]; ++i) {
        r = asprintf(&q,
                     "INSERT INTO ped(BCF_ID, BCF_Sample)"
                     "VALUES (%u, '%s');",
                     i,
                     hdr->samples[i]);

        rc = sqlite3_exec(db, q, NULL, 0, &err_msg);
        if( rc != SQLITE_OK )
            err(EX_SOFTWARE,"SQL error \"%s\" in query \"%s\"", err_msg, q);
    }
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    if (num_ped_fields > 0) {

        fprintf(stderr,
                "Joining values based on BCF_Sample in %s and %s in %s.\n",
                bcf_file_name,
                ped_field_names[col - 1],
                ped_file_name);

        FILE *ped_f = fopen(ped_file_name, "r");
        if (!ped_f)
            err(EX_NOINPUT, "Cannot open file \"%s\"", ped_file_name);

        char *line = NULL;
        size_t len = 0;

        ssize_t read = getline(&line, &len, ped_f); // skip header
        char **ped_values = (char **) malloc(num_ped_fields *sizeof(char *));
        char *tmp_q;

        while ( (read = getline(&line, &len, ped_f)) != -1) {
            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';


            ped_values[0] = strtok(line, "\t");
            for (i = 1; i < num_ped_fields; ++i)
                ped_values[i] = strtok(NULL, "\t");

            r = asprintf(&q, "UPDATE ped SET");

            int first = 0;
            for (i = 0; i < num_ped_fields; ++i) {
                if (first != 0) {
                    r = asprintf(&tmp_q, "%s,", q);
                    free(q);
                    q = tmp_q;
                }

                if (ped_field_types[i] == 1)
                    r = asprintf(&tmp_q,
                                 "%s %s=%s",
                                 q,
                                ped_field_names[i],
                                ped_values[i]);
                else
                     r = asprintf(&tmp_q,
                                 "%s %s='%s'",
                                 q,
                                ped_field_names[i],
                                ped_values[i]);

                free(q);
                q = tmp_q;
                first = 1;
            }

            r = asprintf(&tmp_q,
                         "%s WHERE BCF_Sample == '%s';",
                         q,
                         ped_values[col - 1]);;
            free(q);
            q = tmp_q;

            rc = sqlite3_exec(db, q, NULL, 0, &err_msg);
            if( rc != SQLITE_OK )
                err(EX_SOFTWARE,"SQL error \"%s\" in query \"%s\"", err_msg, q);
        }

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
    char *err_msg;
    int rc = sqlite3_open(ped_db_file, &db);
    if( rc != SQLITE_OK )
        err(EX_NOINPUT,
            "SQL error \"%s\" for database \"%s\"",
            err_msg, ped_db_file);

    char *test_q;
    int r;
    if (strlen(query) == 0)
        r = asprintf(&test_q, "SELECT BCF_ID FROM ped ORDER BY BCF_ID");
    else
        r = asprintf(&test_q,
                     "SELECT BCF_ID FROM ped WHERE %s ORDER BY BCF_ID;",
                     query);

    struct uint32_t_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    rc = sqlite3_exec(db, test_q, uint32_t_ll_callback, &ll, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,"SQL error \"%s\" in query \"%s\"", err_msg, test_q);

    *R = (uint32_t *) malloc(ll.len * sizeof(uint32_t));

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
    char *err_msg;
    int rc = sqlite3_open(ped_db_file, &db);
    if( rc != SQLITE_OK )
        err(EX_NOINPUT,
            "SQL error \"%s\" for database \"%s\"",
            err_msg, ped_db_file);

    char *test_q;
    int r;
    if (strlen(query) == 0)
        r = asprintf(&test_q, "SELECT BCF_ID,%s FROM ped ORDER BY BCF_ID",
                label_id);
    else
        r = asprintf(&test_q,
                     "SELECT BCF_ID,%s FROM ped WHERE %s ORDER BY BCF_ID;",
                     label_id,
                     query);

    struct char_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    rc = sqlite3_exec(db, test_q, char_ll_callback, &ll, &err_msg);
    if( rc != SQLITE_OK )
        err(EX_SOFTWARE,"SQL error \"%s\" in query \"%s\"", err_msg, test_q);

    *R = (char **) malloc(ll.len * sizeof(char *));

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
