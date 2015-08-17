#ifndef __WAH_H__
#define __WAH_H__

#if 0
struct wah_file {
    FILE *file;
    char *file_name;
    uint32_t num_fields, num_records;
    uint32_t word_size;
    uint64_t *record_offsets;
    long header_offset;
    struct gqt_file_header *gqt_header;
};
#endif

struct wah_active_word {
    uint32_t value, nbits;
};

struct wah16_active_word {
    uint16_t value, nbits;
};

struct wah16_ll {
    struct wah16_active_word value;
    struct wah16_ll *next;
};

struct wah_ll {
    struct wah_active_word value;
    struct wah_ll *next;
};

struct wah_run {
    uint32_t *words;
    uint32_t len,
                 word_i,
                 fill, // one word-long version of the fill
                 fill_bit, // if it is a fill, set this bit
                 num_words, // number of words in the run
                 is_fill; //is it a fill run
};

#if 0
/**
 * @brief Open a WAH-encoded (non-bitmap) index and initialize the wah_file
 * data structure.
 *
 * The WAH data structure includes the number of fields, number of records, an
 * array of offsets of the different WAH Memory is allocated for the bitmap
 * offsets within this function and must be freed when its use is complete (
 * free(wf.record_offsets))
 *
 * @param file_name The name of the WAH file
 *
 * @retval WAH data strucutre.
 *
 * Example Usage:
 * @code
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *      struct wah_file wf = init_wah_file(wah_file_name);
 * @endcode
 */
struct wah_file init_wah_file(char *file_name);

void destroy_wah_file(struct wah_file *wf);
#endif

/**
 * @brief   AND two WAH runs
 *
 * @param x WAH run for
 * @param y WAH run for
 * @param O the result of x AND y 
 *
 * @retval  The number of elements in O
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *      uint32_t X[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *      };
 *      uint32_t Y[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *      };
 *      uint32_t *w_X;
 *      int wah_size_X = ints_to_wah(X,5,160,&w_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X);
 *      uint32_t *w_Y;
 *      int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
 *      struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);
 *      uint32_t *Z;
 *      uint32_t Z_len = wah_and(&r_X, &r_Y, &Z);
 * @endcode
 */
uint32_t wah_and(struct wah_run *x,
                     struct wah_run *y,
                     uint32_t **O);


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
 *
 * Example Usage:
 * @code
 *      uint32_t X[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *      };
 *      uint32_t Y[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *      };
 *      uint32_t *w_X;
 *      int wah_size_X = ints_to_wah(X,5,160,&w_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X);
 *      uint32_t *w_Y;
 *      int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
 *      struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);
 *      uint32_t *Z;
 *      uint32_t Z_len = wah_or(&r_X, &r_Y, &Z);
 * @endcode
 */
uint32_t wah_or(struct wah_run *x,
                    struct wah_run *y,
                    uint32_t **O);

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
 *      uint32_t X[5] =
 *              { bin_char_to_int("01000000000000000000000000000001"),
 *                bin_char_to_int("11111111111111111111111111111111"),
 *                bin_char_to_int("11111111111111111111111111111111"),
 *                bin_char_to_int("01000000000101010100000000000000"),
 *                bin_char_to_int("01000000000000000001010101000000")
 *              };
 *
 *      uint32_t *w_X;
 *      uint32_t wah_size_X = ints_to_wah(X,5,160,&w_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X); 
 * @endcode
 */
struct wah_run init_wah_run(uint32_t *words,
                            uint32_t len);

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
 *      uint32_t I[5] = {2147483648,0,0,3,1};
 *      uint32_t *O;
 *      uint32_t wah_size = ints_to_wah(I,5,160,&O);
 *      struct wah_run r = init_wah_run(O, wah_size);
 *      wah_run_decode(&r);
 * @endcode
 */
void wah_run_decode(struct wah_run *r);

/**
 * @brief   Convert an array of 32-bit integers to WAH encoding
 *
 * @param I         An array of 32-bit itergers that enocde a binary string
 * @param I_len     The number of intergers in I
 * @param used_bits The number of used bits (size minus padding)
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 *
 * @ingroup WAH
 *
 * @retval          Size of W
 *
 * Example Usage:
 * @code
 *      uint32_t X[5] =
 *          { bin_char_to_int("01000000000000000000000000000001"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("01000000000101010100000000000000"),
 *            bin_char_to_int("01000000000000000001010101000000")
 *          };
 *      uint32_t *w_X;
 *      uint32_t wah_size_X = ints_to_wah(X,5,160,&w_X);
 * @endcode
 */
uint32_t ints_to_wah(uint32_t *I,
                         int I_len,
                         uint32_t used_bits,
                         uint32_t **W);

/**
 * @brief   Convert an array of 32-bit integers to 16-bit WAH encoding
 *
 * @param I         An array of 32-bit itergers that enocde a binary string
 * @param I_len     The number of intergers in I
 * @param used_bits The number of used bits (size minus padding)
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 *
 * @ingroup WAH
 *
 * @retval          Size of W
 *
 * Example Usage:
 * @code
 *      uint32_t X[5] =
 *          { bin_char_to_int("01000000000000000000000000000001"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("01000000000101010100000000000000"),
 *            bin_char_to_int("01000000000000000001010101000000")
 *          };
 *      uint16_t *w_X;
 *      uint32_t wah_size_X = ints_to_wah(X,5,160,&w_X);
 * @endcode
 */
uint32_t ints_to_wah16(uint32_t *I,
                           int I_len,
                           uint32_t used_bits,
                           uint16_t **W);

/**
 * @brief   A helper function to group ints encoding 32-bits to ones encoding
 *          31-bits for WAH consideration
 *
 * @param I         An array of ints encoding 32-bits
 * @param I_len     Number of elements in I
 * @param used_bits Then number of bits used (size minus padding)
 * @param O         The same bits from I, but encoded in 31-bit groups
 *
 * @ingroup WAH
 *
 * @retval          size of O
 *
 * Example Usage:
 * @code
 *      uint32_t A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      uint32_t *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_31_bits(I,5,160,&O,&l);
 * @endcode
 */
uint32_t map_from_32_bits_to_31_bits(uint32_t *I,
                                         int I_len,
                                         uint32_t used_bits,
                                         uint32_t **O);

/**
 * @brief   A helper function to group ints encoding 32-bits to ones encoding
 *          15-bits for WAH consideration
 *
 * @param I         An array of ints encoding 32-bits
 * @param I_len     Number of elements in I
 * @param used_bits Then number of bits used (size minus padding)
 * @param O         The same bits from I, but encoded in 15-bit groups
 *
 * @ingroup WAH
 *
 * @retval          size of O
 *
 * Example Usage:
 * @code
 *      #include <stdint.h>
 *      uint32_t A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      uint16_t *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_15_bits(I,5,160,&O,&l);
 * @endcode
 */
uint32_t map_from_32_bits_to_15_bits(uint32_t *I,
                                         int I_len,
                                         uint32_t used_bits,
                                         uint16_t **O);

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
 * @brief   Append an active workd to a 16-bit WAH list
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
 *     struct wah16_active_word_ll *A_head = NULL,
 *                               *A_tail = NULL;
 *
 *     struct wah16_active_word a;
 *     a.nbits = 15;
 *     a.value = 19114;
 *     int r = append_active_word(&A_head,&A_tail,a);
 * @endcode
 */
int append_active_16word(struct wah16_ll **A_head,
                         struct wah16_ll **A_tail,
                         struct wah16_active_word a);



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
                     uint32_t fill_size);
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
 *     uint32_t I[5] = {2147483648,0,0,3,1};
 *     uint32_t *WAH;
 *     uint32_t wah_size = ints_to_wah(I,5,160,&WAH);
 *     uint32_t *INTS;
 *     uint32_t ints_size = wah_to_ints(WAH,wah_size,&INTS);
 * @endcode
 */
uint32_t wah_to_ints(uint32_t *W,
                         uint32_t W_len,
                         uint32_t **O);

#if 0
/**
 * @brief Print a WAH encoded (non-bitmap) file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param wf An initilized WAH (non-bitmap) file 
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *                              2:packed int
 *
 *
 * @retval number of records printed 
 */

uint32_t print_wah(struct wah_file wf,
                       uint32_t *record_ids,
                       uint32_t num_r,
                       uint32_t format);
#endif

int append_active_word_b(uint32_t *R,
                         uint32_t R_len,
                         uint32_t value);

uint32_t append_fill_word_b(uint32_t *R,
                            uint32_t R_len,
                            int fill_bit,
                            uint32_t fill_size);

uint32_t wah_or_b(uint32_t *R,
                  uint32_t *X,
                  uint32_t len_X,
                  uint32_t *Y,
                  uint32_t len_Y);

uint32_t uint32_t_to_wah_bitmap(uint32_t *I,
                                uint32_t len_I,
                                uint32_t *bin_range_lo,
                                uint32_t *bin_range_hi,
                                uint32_t len_bin_ranges,
                                uint32_t less_than_bin,
                                uint32_t greater_than_bin,
                                uint32_t num_bins,
                                uint32_t ***wah_bit_maps,
                                uint32_t **wah_bit_map_lens);

uint32_t int_to_wah_bitmap(int *I,
                           uint32_t len_I,
                           int *bin_range_lo,
                           int *bin_range_hi,
                           uint32_t len_bin_ranges,
                           uint32_t less_than_bin,
                           uint32_t greater_than_bin,
                           uint32_t num_bins,
                           uint32_t ***wah_bit_maps,
                           uint32_t **wah_bit_map_lens);

uint32_t float_to_wah_bitmap(float *I,
                             uint32_t len_I,
                             float *bin_range_lo,
                             float *bin_range_hi,
                             uint32_t len_bin_ranges,
                             uint32_t less_than_bin,
                             uint32_t greater_than_bin,
                             uint32_t num_bins,
                             uint32_t ***wah_bit_maps,
                             uint32_t **wah_bit_map_lens);

#endif
