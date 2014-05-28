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

struct wah_ll {
    struct wah_active_word value;
    struct wah_ll *next;
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

/**
 * @brief   Convert an array of 32-bit integers to WAH encoding
 *
 * @param I         An array of 32-bit itergers that enocde a binary string
 * @param I_len     The number of intergers in I
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param W_len     The size of W (set in function)
 *
 * @ingroup WAH
 *
 * @retval          Size of W, will equal W_len
 *
 * Example Usage:
 * @code
 *      unsigned int X[5] =
 *          { bin_char_to_int("01000000000000000000000000000001"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("01000000000101010100000000000000"),
 *            bin_char_to_int("01000000000000000001010101000000")
 *          };
 *      unsigned int *w_X;
 *      int l_X;
 *      int wah_size_X = wah(X,5,&w_X,&l_X);
 * @endcode
 */
int wah(unsigned int *I,
        int I_len,
        unsigned int **W,
        int *W_len);

/**
 * @brief   Add a new bit to an active word
 *
 * @param a
 * @param b
 *
 * @retval  1   Bit could not be added to the active word (it was full, new bit
 *              did not match the fill value)
 * @retval  0   Bit was added to the active word
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *     struct wah_active_word a;
 *     // 1101101 -> 7 bits, int = 109
 *     a.value = 109;
 *     a.nbits = 7;
 *     int r = append_bit_to_active_word(&a, 1);
 * @endcode
 */
int append_bit_to_active_word(struct wah_active_word *a, int b);

/**
 * @brief   Append an active workd to a WAH list
 *
 * @param A_head    Pointer to the head of the list
 * @param A_tail    Pointer to the tail of the list
 * @param a         The active word
 *
 * @retval  1   A new word was added to the list
 * @retval  0   A new word was not added
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *     struct wah_active_word_ll *A_head = NULL,
 *                               *A_tail = NULL;
 *
 *     struct wah_active_word a;
 *     a.nbits = 33;
 *     a.value = 2147483681;
 *     int r = append_active_word(&A_head,&A_tail,a);
 * @endcode
 */
int append_active_word(struct wah_ll **A_head,
                       struct wah_ll **A_tail,
                       struct wah_active_word a);
/**
 * @brief   Append a fill word to a WAH list
 *
 * @param A_head    Pointer to the head of the list
 * @param A_tail    Pointer to the tail of the list
 * @param fill_bit  Value of the fill (1 or 0)
 * @param fill_size Number of words in the fill
 *
 * @retval  1   A new word was added to the list
 * @retval  0   A new word was not added
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *     struct wah_active_word_ll *A_head = NULL,
 *                               *A_tail = NULL;
 *     struct wah_active_word a;
 *     int r = append_fill_word(&A_head,&A_tail,1,1);
 * @endcode
 */
int append_fill_word(struct wah_ll **A_head,
                     struct wah_ll **A_tail,
                     int fill_bit,
                     unsigned int fill_size);
/**
 * @brief   Initialize the values of WAH run for subsequent binary op
 *
 * @param words An array of WAH encoded values
 * @param len   Number of elements in words
 *
 * @retval  run associated with the values in words
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *      unsigned int X[5] =
 *              { bin_char_to_int("01000000000000000000000000000001"),
 *                bin_char_to_int("11111111111111111111111111111111"),
 *                bin_char_to_int("11111111111111111111111111111111"),
 *                bin_char_to_int("01000000000101010100000000000000"),
 *                bin_char_to_int("01000000000000000001010101000000")
 *              };
 *
 *      unsigned int *w_X;
 *      int l_X;
 *      int wah_size_X = wah(X,5,&w_X,&l_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X); 
 * @endcode
 */
struct wah_run init_wah_run(unsigned int *words,
                            unsigned int len);

/**
 * @brief   Decoded the word-ith element
 *
 * Take the word_i-th element in the run. If that word is a fill, set fill
 * to all 1s or 0s based on the fill bit, num_words to the number of words that
 * the fill covers, and is_fill to 1.  If the word-i-th element is a litteral,
 * set is_fill to 0
 *
 * @param r pointer to a WAH run
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *      unsigned int I[5] = {2147483648,0,0,3,1};
 *      unsigned int *O;
 *      int l;
 *      int wah_size = wah(I,5,&O,&l);
 *      struct wah_run r = init_wah_run(O, wah_size);
 *      wah_run_decode(&r);
 * @endcode
 */
void wah_run_decode(struct wah_run *r);

/**
 * @brief   OR two WAH runs
 *
 * @ingroup WAH
 */
void wah_or(struct wah_run *x, struct wah_run *y);


/**
 * @brief   A helper function for testing
 *
 * @param I         An array of ints encoding 32-bits
 * @param I_len     Number of elements in I
 * @param O         The same bits from I, but encoded in 31-bit groups
 * @param O_len     Number of elements in O
 *
 * @ingroup TESTING
 *
 * @retval          size of O, will equal O_len
 *
 * Example Usage:
 * @code
 *      unsigned int A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      unsigned int *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_31_bits(I,5,&O,&l);
 * @endcode
 */
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
