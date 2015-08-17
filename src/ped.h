#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include "genotq.h"

#ifndef __PED_H__
#define __PED_H__

struct uint32_t_ll_node {
    uint32_t v;
    struct uint32_t_ll_node *next;
};

struct uint32_t_ll {
    struct uint32_t_ll_node *head, *tail;
    uint32_t len;
};


struct char_ll_node {
    char *v;
    struct char_ll_node *next;
};

struct char_ll {
    struct char_ll_node *head, *tail;
    uint32_t len;
};


static int uint32_t_ll_callback(void *ll_p,
                                int argc,
                                char **argv,
                                char **col_name);

uint32_t convert_file_by_name_ped_to_db(char *bcf_file_name,
                                        char *ped_file_name,
                                        uint32_t col,
                                        char *db_name);

uint32_t resolve_ind_query(uint32_t **R, char *query, char *ped_db_file);

uint32_t resolve_label_query(char ***R,
                             char *label_id,
                             char *query,
                             char *ped_db_file);
#endif
