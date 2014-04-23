#ifndef __GENOTQ_H__
#define __GENOTQ_H__

struct ubin_file {
    FILE *file;
    unsigned int num_fields, num_records;
};

struct ubin_file init_ubin_file(char *file_name);

void init_int_genotype_reader(char *file_name, int num_gt);

void destroy_int_genotype_reader();

int get_next_int_genotype(int *line_num, int *gt_num, int *gt);

int get_ubin_genotypes(struct ubin_file u_file,
                       int num_gts,
                       int *record_ids,
                       int *field_ids,
                       int *gts);

int or_ubin_records(struct ubin_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r);

int or_ubin_fields(struct ubin_file u_file, 
                    int num_f,
                    int *field_ids,
                    unsigned int **r);

void parse_cmd_line_int_csv(int *I,
                            int num_I,
                            char *cmd_line_arg);
#endif
