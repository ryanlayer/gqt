#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include "genotq.h"

//{{{ static int callback(void *ll_p,
static int callback(void *ll_p,
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

//{{{ uint32_t convert_file_by_name_ped_to_db(char *ped_file_name, char
uint32_t convert_file_by_name_ped_to_db(char *ped_file_name, char *db_name)
{
    // remove the db if it is already there, replace this later
    struct stat buffer;   
    int r = stat(db_name, &buffer);

    if (r == 0)
        remove(db_name);


    FILE *ped_f = fopen(ped_file_name, "r");
    if (ped_f == NULL)
        exit(EXIT_FAILURE);


    char *line = NULL, *tmp_line;
    size_t len = 0;

    ssize_t read = getline(&line, &len, ped_f);
    if (line[strlen(line) - 1] == '\n')
        line[strlen(line) - 1] = '\0';

    char *word;
    uint32_t num_fields = 0;

    tmp_line = (char *) malloc(strlen(line) * sizeof(char));
    strcpy(tmp_line, line);

    word = strtok(tmp_line, "\t");
    while (word != NULL) {
        num_fields += 1;
        word = strtok(NULL, "\t");
    }

    num_fields += 1; // one more field for the ind_id

    char **field_name = (char **) malloc (num_fields * sizeof(char *));

    strcpy(tmp_line, line);

    word = strtok(tmp_line, "\t");
    uint32_t i = 0, j;
    while (word != NULL) {
        char *name = (char *) malloc(strlen(word)*sizeof(char));
        strcpy(name, word);
        for (j = 0; j < strlen(name); ++j) {
            if (name[j] == ' ')
                name[j] = '_';
        }
        field_name[i] = name;
        i += 1;
        word = strtok(NULL, "\t");
    }

    char *ind_id_name = "Ind_ID";
    field_name[i] = ind_id_name;

    int *field_type = (int *) malloc(num_fields * sizeof(int));
    for (i = 0; i < num_fields; ++i) 
        field_type[i] = 1;

    while ( (read = getline(&line, &len, ped_f)) != -1) {
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        uint32_t j;
        word = strtok(line, "\t");
        for (i = 0; i < num_fields - 1; ++i) { // skip added last ind_id field
            for (j = 0; j < strlen(word); ++j) 
                field_type[i] &= isnumber((int)word[j]);
            word = strtok(NULL, "\t");
        }
    }
    fclose(ped_f);

    char *q_create_table, *q_create_table_tmp;
    r = asprintf(&q_create_table, "CREATE TABLE ped(");
    for (i = 0; i < num_fields; ++i) {
        if (i != 0) {
            r = asprintf(&q_create_table_tmp, "%s ", q_create_table);
            free(q_create_table);
            q_create_table = q_create_table_tmp;
        }

        if (field_type[i] == 1)
            r = asprintf(&q_create_table_tmp, 
                         "%s %s INTEGER",
                         q_create_table,
                         field_name[i]);
                         
        else
            r = asprintf(&q_create_table_tmp,
                         "%s %s TEXT",
                         q_create_table,
                         field_name[i]);
        free(q_create_table);
        q_create_table = q_create_table_tmp;

        if (i < (num_fields - 1)) {
            r = asprintf(&q_create_table_tmp, "%s,", q_create_table);
            free(q_create_table);
            q_create_table = q_create_table_tmp;
        }
    }

    r = asprintf(&q_create_table_tmp, "%s);", q_create_table);
    free(q_create_table);
    q_create_table = q_create_table_tmp;


    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(db_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
    }

    rc = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
    }

    char *i_pre_text, *i_pre_tmp;
    r = asprintf(&i_pre_text, "INSERT INTO ped(");
    for (i = 0; i < num_fields; ++i) {
        if (i != 0) {
            r = asprintf(&i_pre_tmp, "%s,", i_pre_text);
            free(i_pre_text);
            i_pre_text = i_pre_tmp;
        }
        r = asprintf(&i_pre_tmp, "%s%s", i_pre_text, field_name[i]);
        free(i_pre_text);
        i_pre_text = i_pre_tmp;
    }
    r = asprintf(&i_pre_tmp, "%s)", i_pre_text);
    free(i_pre_text);
    i_pre_text = i_pre_tmp;


    ped_f = fopen(ped_file_name, "r");
    if (ped_f == NULL)
        exit(EXIT_FAILURE);

    // get rid of the header line
    read = getline(&line, &len, ped_f);

    char *q_insert, *q_insert_tmp;
    uint32_t ind_id = 0;
    while ( (read = getline(&line, &len, ped_f)) != -1) {
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        r = asprintf(&q_insert, "%s VALUES(", i_pre_text);

        word = strtok(line, "\t");
        for (i = 0; i < num_fields - 1; ++i) { // skip last ind_id field
            if ( i != 0 ) {
                r = asprintf(&q_insert_tmp, "%s,", q_insert);
                free(q_insert);
                q_insert = q_insert_tmp;
            }
            if (field_type[i] == 1)
                r = asprintf(&q_insert_tmp, "%s %s", q_insert, word);
            else
                r = asprintf(&q_insert_tmp, "%s '%s'", q_insert, word);
            free(q_insert);
            q_insert = q_insert_tmp;
            word = strtok(NULL, "\t");
        }
        r = asprintf(&q_insert_tmp, "%s,%d);", q_insert, ind_id);
        free(q_insert);
        q_insert = q_insert_tmp;

        rc = sqlite3_exec(db, q_insert, NULL, 0, &err_msg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", err_msg);
            fprintf(stderr, "%s\n", q_insert);
            sqlite3_free(err_msg);
        }

        ind_id += 1;
        free(q_insert);
    }
    fclose(ped_f);
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
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
    }

    char *test_q;
    int r = asprintf(&test_q, "SELECT Ind_ID FROM ped WHERE %s;", query);

    struct uint32_t_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    rc = sqlite3_exec(db, test_q, callback, &ll, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
    }

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
