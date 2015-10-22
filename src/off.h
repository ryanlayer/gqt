#ifndef __OFF_H__
#define __OFF_H__

#include <stdint.h>
#include <stdio.h>

struct off_file {
    FILE *file;
    char *file_name;
    struct gqt_file_header *gqt_header;
    uint64_t *offsets;
};

struct off_file *new_off_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples);

void add_to_off_file(struct off_file *b, uint64_t next_off);

struct off_file *open_off_file(char *file_name);

void destroy_off_file(struct off_file *b);

#endif
