#ifndef __GENOTQ_H__
#define __GENOTQ_H__

struct plt_file {
    FILE *file;
    unsigned int num_fields, num_records;
    long header_offset;
};

struct ubin_file {
    FILE *file;
    unsigned int num_fields, num_records;
    long header_offset;
};

struct uint_ll {
    unsigned int value;
    struct uint_ll *next;
};

struct wah_active_word {
    unsigned int value, nbits;
};

struct wah_active_word_ll {
    struct wah_active_word value;
    struct wah_active_word_ll *next;
};

struct wah_run {
    unsigned int *words;
    unsigned int len,
                 word_i,
                 fill, // one word-long version of the fill
                 num_words, // number of words in the run
                 is_fill; //is it a fill run
};

unsigned int bin_char_to_int(char *bin);

struct wah_run init_wah_run(unsigned int *words,
                            unsigned int len);

void wah_run_decode(struct wah_run *r);

void wah_or(struct wah_run *x, struct wah_run *y);

struct plt_file init_plt_file(char *file_name);

int or_records_plt(struct plt_file pf,
                   int *record_ids,
                   int num_r,
                   int *G);

int or_fields_plt(struct plt_file pf,
                  int *field_ids,
                  int num_f,
                  int *G);

unsigned int pack_2_bit_ints(int *ints, int num_ints);

int *unpack_2_bit_ints(unsigned int packed_ints);


int plt_by_name_to_ubin(char *in_file_name, char *out_file_name);

int plt_to_ubin(struct plt_file pf, char *out_file_name);

int plt_line_to_packed_ints(char *line,
                            int num_fields, 
                            unsigned int **packed_ints,
                            int *len);

int or_records_ubin(struct ubin_file uf, 
                    int *record_ids,
                    int num_r,
                    unsigned int **G);

int or_fields_ubin(struct ubin_file uf, 
                   int *fields_ids,
                   int num_f,
                   unsigned int **G);


unsigned int rle(unsigned int *I, int I_len, unsigned int **O);

int wah(unsigned int *I,
        int I_len,
        unsigned int **O,
        int *O_len);

int append_bit_to_active_word(struct wah_active_word *a, int b);

int append_active_word(struct wah_active_word_ll **A_head,
                       struct wah_active_word_ll **A_tail,
                       struct wah_active_word a);

int map_from_32_bits_to_31_bits(unsigned int *I,
                                int I_len,
                                unsigned int **O,
                                int *O_len);

////////////////////////////////////////////////////////////////////////


struct uint_file {
    FILE *file;
    unsigned int num_fields, num_records;
    int line_len;
};




struct uint_file init_uint_file(char *file_name, 
                                int num_records,
                                int num_fields);

int or_uint_records(struct uint_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r);

int or_uint_fields(struct uint_file u_file, 
                    int num_f,
                    int *field_ids,
                    unsigned int **r);

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
