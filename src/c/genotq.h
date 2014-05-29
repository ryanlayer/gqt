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
                 fill_bit, // if it is a fill, set this bit
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


int convert_file_plt_by_name_to_ubin(char *in_file_name, char *out_file_name);

int convert_file_plt_to_ubin(struct plt_file pf, char *out_file_name);

unsigned int plt_line_to_packed_ints(char *line,
                                     int num_fields, 
                                     unsigned int **packed_ints);

int or_records_ubin(struct ubin_file uf, 
                    int *record_ids,
                    int num_r,
                    unsigned int **G);

int or_fields_ubin(struct ubin_file uf, 
                   int *fields_ids,
                   int num_f,
                   unsigned int **G);


/**
 * @brief   Compress an array of integers encoded binary digits using
 *          run-length encoding.
 *
 * @param I     Array of integers
 * @param I_len Number of elements in I
 * @param O     Array containing the run-length encoding of I
 *
 * @retval number of elements in O
 */
unsigned int ints_to_rle(unsigned int *I, int I_len, unsigned int **O);

/**
 * @brief   Convert an array of 32-bit integers to WAH encoding
 *
 * @param I         An array of 32-bit itergers that enocde a binary string
 * @param I_len     The number of intergers in I
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 *
 * @ingroup WAH
 *
 * @retval          Size of W
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
 *      unsigned int wah_size_X = ints_to_wah(X,5,&w_X);
 * @endcode
 */
unsigned int ints_to_wah(unsigned int *I,
                         int I_len,
                         unsigned int **W);

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
 *      unsigned int wah_size_X = ints_to_wah(X,5,&w_X);
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
 *      unsigned int wah_size = ints_to_wah(I,5,&O);
 *      struct wah_run r = init_wah_run(O, wah_size);
 *      wah_run_decode(&r);
 * @endcode
 */
void wah_run_decode(struct wah_run *r);

/**
 * @brief   OR two WAH runs
 *
 * @param x WAH run for
 * @param y WAH run for
 * @param O the result of x OR y 
 *
 * @retval  The number of elements in O
 *
 * @ingroup WAH
 */
unsigned int wah_or(struct wah_run *x,
                    struct wah_run *y,
                    unsigned int **O);


/**
 * @brief   A helper function to group ints encoding 32-bits to ones encoding
 *          31-bits for WAH consideration
 *
 * @param I         An array of ints encoding 32-bits
 * @param I_len     Number of elements in I
 * @param O         The same bits from I, but encoded in 31-bit groups
 *
 * @ingroup WAH
 *
 * @retval          size of O
 *
 * Example Usage:
 * @code
 *      unsigned int A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      unsigned int *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_31_bits(I,5,&O,&l);
 * @endcode
 */
unsigned int map_from_32_bits_to_31_bits(unsigned int *I,
                                         int I_len,
                                         unsigned int **O);

/**
 * @brief Convert WAH encoding to uncompressed binary ints.
 *
 * Given that padding must be used to deal with the  31-bit chunks that WAH
 * considers and the 32-bit chunks of an int, it is common to have more bits in
 * the resulting int than were in the original.
 *
 * @param W     A WAH-encoded array
 * @param W_len The number of elements in W
 * @param O     The uncompressed binary version of W
 *
 * @retval size of O
 *
 * Example Usage:
 * @code
 *     unsigned int I[5] = {2147483648,0,0,3,1};
 *     unsigned int *WAH;
 *     unsigned int wah_size = ints_to_wah(I,5,&WAH);
 *     unsigned int *INTS;
 *     unsigned int ints_size = wah_to_ints(WAH,wah_size,&INTS);
 * @endcode
 */
unsigned int wah_to_ints(unsigned int *W,
                         unsigned int W_len,
                         unsigned int **O);

/**
 * @brief Convert an array of uncompressed binary values to a bitmap index of
 *        the values
 *
 * @param U Uncompressed 2-bit binary values encoded in an array of integers
 * @param U_len Unumber of integers in U
 * @param B Resulting bitmap index where the index for value 0 starts at
 *          position 0 and the index for value 1 starts at position N, value 2
 *          at N*2, etc. Where U encodeds N (U_len*16) 2-bit values
 * 
 * @retval size total number of elements in B, where there bitmap index for each
 *         value uses 1/4 of the size
 *
 * Example Usage:
 * @code
 *      unsigned int U1[4] = {
 *              bin_char_to_int("00011011000110110001101100011011"),
 *              bin_char_to_int("00000101010110101111000001010110"),
 *              bin_char_to_int("00000000000000000000000000000000"),
 *              bin_char_to_int("11111111111111111111111111111111")
 *      };
 *      unsigned int *B;
 *      unsigned int B_len = ubin_to_bitmap(U1,4,u,&B);
 * @endcode
 */
unsigned int  ubin_to_bitmap(unsigned int *U,
                             unsigned int U_len,
                             unsigned int **B);


/**
 * @brief Convert an array of uncompressed two-bit binary values to WAH
 *        encoding of the bitmap index representation of the two-bit binary
 *
 *  In this scheme, there are 4 uniq values in the uncompressed binary
 *  (0,1,2,3), so the WAH encoding will create 4 different sets, one for each
 *  bitmap index
 *
 * @param U         Uncompressed 2-bit binary values encoded in an array of
 *                  integers
 * @param U_len     Unumber of integers in U
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param offsets   A list of 4 offsets (indicies within W) to the starting
 *                  locations of each bitmap index.  This size of each index
 *                  can be found by takeing the differnet of the offsets and
 *                  the total size (returned by fuction)
 *
 * @retval number of ints in W
 *
 * Example Usage:
 * @code
 *    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";
 *
 *    unsigned int *ubin;
 *    unsigned int ubin_len = plt_line_to_packed_ints(plt, 192, &ubin);
 *
 *    TEST_ASSERT_EQUAL(12, ubin_len);
 *
 *    unsigned int *wah;
 *    unsigned int *wah_offsets;
 *    unsigned int wah_len = ubin_to_bitmap_wah(ubin,
 *                                              ubin_len,
 *                                              &wah,
 *                                              &wah_offsets);
 * @endcode
 */
unsigned int ubin_to_bitmap_wah(unsigned int *U,
                                unsigned int U_len,
                                unsigned int **W,
                                unsigned int **offsets);
/**
 * @brief Convert an string of plain text values over the alphabet [0,1,2,3]
 *        to WAH encoding of the bitmap index representation of the two-bit
 *        binary
 *
 *  In this scheme, there are 4 uniq values in plain text (0,1,2,3), so the WAH
 *  encoding will create 4 different sets, one for each bitmap index
 *
 * @param P         Space-seperated string of 0,1,2, or 3
 * @param P_len     Number of digits in P
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param offsets   A list of 4 offsets (indicies within W) to the starting
 *                  locations of each bitmap index.  This size of each index
 *                  can be found by takeing the differnet of the offsets and
 *                  the total size (returned by fuction)
 *
 * @retval number of ints in W
 *
 * Example Usage:
 * @code
 *
 * @endcode
 */

unsigned int plt_to_bitmap_wah(char *plt,
                               unsigned int plt_len,
                               unsigned int **W,
                               unsigned int **offsets);


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
