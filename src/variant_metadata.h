#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>


#ifndef __VARIANT_METADATA_H__
#define __VARIANT_METADATA_H__

struct variant_md_info {
    int rowid;
    uint64_t offset;
    uint32_t bins;
    int type;
    char *source_file;
};

struct bin_info {
    uint32_t bin, is_hi, is_lo;
    float lo, hi;
    int was_set;
};

int register_variant_metadata_index(char *bcf_file_name,
                                    char *db_file_name,
                                    char *field_name,
                                    uint64_t start_offset,
                                    int type,
                                    uint32_t num_bins);

uint32_t add_variant_metadata_float_bins(char *db_file_name,
                                         int rowid,
                                         float *float_bin_range_lo,
                                         float *float_bin_range_hi,
                                         int actual_num_bins,
                                         int less_than_bin,
                                         int greater_than_bin);
int count_callback(void *count,
                          int argc,
                          char **argv,
                          char **col_name);

/*
 * Use an op flag
 * 1: equal
 * 2: less_than
 * 4: less_than_equal
 * 8: greater_than
 * 16: greater_than_equal
 * 32: not_equal
 *
 * if 1 or 32 are set the lo value is used
 * if 2 or 4 are set lo value is used
 * if 8 or 16 are set hi value is used
 */
int get_variant_metadata_query_bins(char *variant_db_name,
                                    char *field_name,
                                    int rowid,
                                    uint32_t num_bins,
                                    int type,
                                    int op_flag,
                                    float lo,
                                    float hi,
                                    uint32_t **bins);

int get_variant_metadata_bin_info(char *variant_db_name,
                                  char *field_name,
                                  uint64_t *offset,
                                  uint32_t *num_bins,
                                  int *type,
                                  char **source_file);

int get_single_int_callback(void *val,
                            int argc,
                            char **argv,
                            char **col_name);

#endif
