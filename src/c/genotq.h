#ifndef __GENOTQ_H__
#define __GENOTQ_H__

#include <stdint.h>


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

struct wah_file {
    FILE *file;
    unsigned int num_fields, num_records;
    unsigned int *record_offsets;
    long header_offset;
};

struct uint_ll {
    unsigned int value;
    struct uint_ll *next;
};

struct wah_active_word {
    unsigned int value, nbits;
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

int *unpack_1_bit_ints(unsigned int packed_ints);

/**
 * @brief Convert a VCF file to a by-variant plain text file
 *
 * @param in_file_name VCF file name
 * @param num_fields Number of fields in the VCF
 * @param num_records Number of records in the VCF
 * @param out_file_name plain text file name
 *
 * @retval 0 if things worked
 * @retval 1 otherwise
 *
 * Example Usage:
 * @code
 * @endcode
 */
int convert_file_by_name_vcf_to_plt(char *in_file_name,
                                    unsigned int num_fields,
                                    unsigned int num_records,
                                    char *out_file_name);

/**
 * @brief Convert a plain text file to an uncompressed binary file
 *
 * @param in_file_name Plain text file name
 * @param out_file_name Uncompressed binary file name
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.ind.txt";
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *
 *     int r = convert_file_by_name_plt_to_ubin(plt_file_name,ubin_file_name);
 * @endcode
 */
int convert_file_by_name_plt_to_ubin(char *in_file_name, char *out_file_name);

/**
 * @brief Convert a plain text file to an uncompressed binary file
 *
 * @param pf initialized plt file data structure
 * @param out_file_name name of the uncompressed binary file
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.ind.txt";
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *
 *     struct plt_file pf = init_plt_file(plt_file_name);
 *     int r = convert_file_plt_to_ubin(pf, ubin_file_name);
 *     fclose(pf.file);
 * @endcode
 */
int convert_file_plt_to_ubin(struct plt_file pf, char *out_file_name);

unsigned int plt_line_to_packed_ints(char *line,
                                     unsigned int num_fields, 
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
 * @brief Get a pointer to the uncompressed binary encoded record
 *
 * @param uf The ubin data structure
 * @param record_id The index of the record of interest
 * @param ubin_record A pointer to the uncompressed binary record
 *
 * @retval number of integers in the record
 *
 * Example Usage:
 * @code
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *     struct ubin_file uf = init_ubin_file(ubin_file_name);
 *
 *     unsigned int *ints, num_ints;
 *     unsigned int record_id = 0;
 *     num_ints = get_ubin_record(uf, record_id, &ints);
 * @endcode
 */
unsigned int get_ubin_record(struct ubin_file uf,
                             unsigned int record_id,
                             unsigned int **ubin_record);

/**
 * @brief   Compress an array of integers encoded binary digits using
 *          run-length encoding.
 *
 * @param I     Array of integers
 * @param I_len Number of elements in I
 * @param O     Array containing the run-length encoding of I
 *
 * @retval number of elements in O
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int ints_to_rle(unsigned int *I, int I_len, unsigned int **O);

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
 *      unsigned int X[5] =
 *          { bin_char_to_int("01000000000000000000000000000001"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("01000000000101010100000000000000"),
 *            bin_char_to_int("01000000000000000001010101000000")
 *          };
 *      unsigned int *w_X;
 *      unsigned int wah_size_X = ints_to_wah(X,5,160,&w_X);
 * @endcode
 */
unsigned int ints_to_wah(unsigned int *I,
                         int I_len,
                         unsigned int used_bits,
                         unsigned int **W);

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
 *      unsigned int X[5] =
 *          { bin_char_to_int("01000000000000000000000000000001"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("11111111111111111111111111111111"),
 *            bin_char_to_int("01000000000101010100000000000000"),
 *            bin_char_to_int("01000000000000000001010101000000")
 *          };
 *      uint16_t *w_X;
 *      unsigned int wah_size_X = ints_to_wah(X,5,160,&w_X);
 * @endcode
 */
unsigned int ints_to_wah16(unsigned int *I,
                           int I_len,
                           unsigned int used_bits,
                           uint16_t **W);

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
 *      unsigned int wah_size_X = ints_to_wah(X,5,160,&w_X);
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
 *      unsigned int wah_size = ints_to_wah(I,5,160,&O);
 *      struct wah_run r = init_wah_run(O, wah_size);
 *      wah_run_decode(&r);
 * @endcode
 */
void wah_run_decode(struct wah_run *r);

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
 *      unsigned int X[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *      };
 *      unsigned int Y[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *      };
 *      unsigned int *w_X;
 *      int wah_size_X = ints_to_wah(X,5,160,&w_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X);
 *      unsigned int *w_Y;
 *      int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
 *      struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);
 *      unsigned int *Z;
 *      unsigned int Z_len = wah_and(&r_X, &r_Y, &Z);
 * @endcode
 */
unsigned int wah_and(struct wah_run *x,
                     struct wah_run *y,
                     unsigned int **O);


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
 *      unsigned int X[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *      };
 *      unsigned int Y[5] = {
 *          bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *      };
 *      unsigned int *w_X;
 *      int wah_size_X = ints_to_wah(X,5,160,&w_X);
 *      struct wah_run r_X = init_wah_run(w_X, wah_size_X);
 *      unsigned int *w_Y;
 *      int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
 *      struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);
 *      unsigned int *Z;
 *      unsigned int Z_len = wah_or(&r_X, &r_Y, &Z);
 * @endcode
 */
unsigned int wah_or(struct wah_run *x,
                    struct wah_run *y,
                    unsigned int **O);

/**
 * @brief   in-place OR two WAH runs
 *
 * The assumption here is that r_wah is pre-allocated and  contains only
 * literals, no fills.  Wah can be any combination of litterals and fills.
 * The result is stored back into r_wah.
 *
 * @param r_wah a wah-encoded bitmap of all literals
 * @param r_wah_size the number of words in r_wah
 * @param wah a wah-encoded bitmap with a cobmo of literals and fills
 * @param wah_size the numberof words in wah
 *
 * @retval  The number of elements in r_wah
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *    unsigned int X[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *        };
 *
 *    unsigned int Y[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *        };
 *    unsigned int R[6] = {0,0,0,0,0,0};
 *
 *    unsigned int r = wah_in_place_or(R, 6, w_Y, wah_size_Y);
 *    r = wah_in_place_or(R, 6, w_X, wah_size_X);
 *    unsigned int *ints_ip;
 *    unsigned int ints_ip_len = wah_to_ints(R,r,&ints_ip);
 * @endcode
 */
unsigned int  wah_in_place_or(unsigned int *r_wah,
                              unsigned int r_wah_size,
                              unsigned int *wah,
                              unsigned int wah_size);

/**
 * @brief   in-place AND two WAH runs
 *
 * The assumption here is that r_wah is pre-allocated and  contains only
 * literals, no fills.  Wah can be any combination of litterals and fills.
 * The result is stored back into r_wah.
 *
 * @param r_wah a wah-encoded bitmap of all literals
 * @param r_wah_size the number of words in r_wah
 * @param wah a wah-encoded bitmap with a cobmo of literals and fills
 * @param wah_size the numberof words in wah
 *
 * @retval  The number of elements in r_wah
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 *    unsigned int X[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *        };
 *
 *    unsigned int Y[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *        };
 *    unsigned int R[6] = {0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff};
 *
 *    unsigned int r = wah_in_place_and(R, 6, w_Y, wah_size_Y);
 *    r = wah_in_place_and(R, 6, w_X, wah_size_X);
 *    unsigned int *ints_ip;
 *    unsigned int ints_ip_len = wah_to_ints(R,r,&ints_ip);
 * @endcode
 */
unsigned int  wah_in_place_and(unsigned int *r_wah,
                               unsigned int r_wah_size,
                               unsigned int *wah,
                               unsigned int wah_size);

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
 *      unsigned int A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      unsigned int *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_31_bits(I,5,160,&O,&l);
 * @endcode
 */
unsigned int map_from_32_bits_to_31_bits(unsigned int *I,
                                         int I_len,
                                         unsigned int used_bits,
                                         unsigned int **O);

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
 *      unsigned int A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
 *      uint16_t *O;
 *      int l;
 *      int num_31_groups = map_from_32_bits_to_15_bits(I,5,160,&O,&l);
 * @endcode
 */
unsigned int map_from_32_bits_to_15_bits(unsigned int *I,
                                         int I_len,
                                         unsigned int used_bits,
                                         uint16_t **O);

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
 *     unsigned int wah_size = ints_to_wah(I,5,160,&WAH);
 *     unsigned int *INTS;
 *     unsigned int ints_size = wah_to_ints(WAH,wah_size,&INTS);
 * @endcode
 */
unsigned int wah_to_ints(unsigned int *W,
                         unsigned int W_len,
                         unsigned int **O);

/**
 * @brief Convert an array of uncompressed binary values to a bitmap index of
 *        the values.
 *
 * @param U Uncompressed 2-bit binary values encoded in an array of integers
 * @param U_len Unumber of integers in U
 * @param used_bits Number of bits in U that are used
 * @param B Resulting bitmap index where the index for value 0 starts at
 *          position 0 and the index for value 1 starts at position N, value 2
 *          at N*2, etc. Where U encodeds N (U_len*16) 2-bit values
 * 
 * @retval size total number of elements in B, where the bitmap index for each
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
 *      unsigned int B_len = ubin_to_bitmap(U1,4,128,u,&B);
 * @endcode
 */
unsigned int  ubin_to_bitmap(unsigned int *U,
                             unsigned int U_len,
                             unsigned int used_bits,
                             unsigned int **B);

/**
 * @brief Convert an array of uncompressed binary values to a bitmap index of
 *        the values using 16-bit ints.
 *
 * @param U Uncompressed 2-bit binary values encoded in an array of integers
 * @param U_len Unumber of integers in U
 * @param used_bits Number of bits in U that are used
 * @param B Resulting bitmap index where the index for value 0 starts at
 *          position 0 and the index for value 1 starts at position N, value 2
 *          at N*2, etc. Where U encodeds N (U_len*16) 2-bit values
 * 
 * @retval size total number of elements in B, where the bitmap index for each
 *         value uses 1/4 of the size
 *
 * Example Usage:
 * @code
 *      #include <stdint.h>
 *      unsigned int U1[4] = {
 *              bin_char_to_int("00011011000110110001101100011011"),
 *              bin_char_to_int("00000101010110101111000001010110"),
 *              bin_char_to_int("00000000000000000000000000000000"),
 *              bin_char_to_int("11111111111111111111111111111111")
 *      };
 *      uint16_t *B;
 *      unsigned int B_len = ubin_to_bitmap(U1,4,128,u,&B);
 * @endcode
 */

unsigned int ubin_to_bitmap_wah16(unsigned int *U,
                                  unsigned int U_len,
                                  unsigned int num_fields,
                                  uint16_t **W,
                                  unsigned int **wah_sizes);


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
 * @param num_fields    The number of fields
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param wah_sizes A list of the 4 sizes 
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
 *    unsigned int *wah_sizes;
 *    unsigned int wah_len = ubin_to_bitmap_wah(ubin,
 *                                              ubin_len,
 *                                              192,
 *                                              &wah,
 *                                              &wah_sizes);
 * @endcode
 */
unsigned int ubin_to_bitmap_wah(unsigned int *U,
                                unsigned int U_len,
                                unsigned int num_fields,
                                unsigned int **W,
                                unsigned int **wah_sizes);
/**
 * @brief Convert an string of plain text values over the alphabet [0,1,2,3]
 *        to WAH encoding of the bitmap index representation of the two-bit
 *        binary
 *
 *  In this scheme, there are 4 uniq values in plain text (0,1,2,3), so the WAH
 *  encoding will create 4 different sets, one for each bitmap index
 *
 * @param plt         Space-seperated string of 0,1,2, or 3
 * @param plt_len     Number of digits in P
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param wah_sizes   A list the 4 sizes
 *
 * @retval number of ints in W
 *
 * Example Usage:
 * @code
 *      char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";
 *
 *      unsigned int *wah;
 *      unsigned int *wah_sizes;
 *      unsigned int wah_len = plt_to_bitmap_wah(plt,
 *                                               192,
 *                                               &wah,
 *                                               &wah_sizes);
 * @endcode
 */
unsigned int plt_to_bitmap_wah(char *plt,
                               unsigned int plt_len,
                               unsigned int **W,
                               unsigned int **wah_sizes);


/**
 * @brief convert an uncompressed binary file to a WAH encoded bitmap index
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 4 * 32-bit)
 * 1st record 0 bitmap
 * 1st record 1 bitmap
 * 1st record 2 bitmap
 * 1st record 3 bitmap
 * 2nd record 0 bitmap
 * 2nd record 1 bitmap
 * 2nd record 2 bitmap
 * 2nd record 3 bitmap
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out WAH encoded bitmap index file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wahbm(ubin_file_name, wah_file_name);
 * @endcode
 */
unsigned int convert_file_by_name_ubin_to_wahbm(char *ubin_in, char *wah_out);

/**
 * @brief convert an uncompressed binary file to a 16-bit WAH encoded bitmap
 * index
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 4 * 16-bit)
 * 1st record 0 bitmap
 * 1st record 1 bitmap
 * 1st record 2 bitmap
 * 1st record 3 bitmap
 * 2nd record 0 bitmap
 * 2nd record 1 bitmap
 * 2nd record 2 bitmap
 * 2nd record 3 bitmap
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out 16-bit WAH encoded bitmap index file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah16";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wahbm(ubin_file_name, wah_file_name);
 * @endcode
 */
unsigned int convert_file_by_name_ubin_to_wahbm16(char *ubin_in,
                                                  char *wah_out);
/**
 * @brief Open a WAH-encoded bitmap index and initialize the wah_file data
 * structure.
 *
 * The WAH data structure includes the number of fields, number of records, an
 * array of offsets of the different WAH bitmaps, and an offset of the first
 * bitmap.  Memory is allocated for the bitmap offsets within this function and
 * must be freed when its use is complete ( free(wf.record_offsets))
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
struct wah_file init_wahbm_file(char *file_name);

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

/**
 * @brief Get a pointer to the bitmap of a particular WAH-encoded bitmap record
 *
 * @param wf The WAH file data structure
 * @param wah_record The record ID
 * @param bitmap The bitmap ID (0,1,2, or 3)
 * @param wah_bitmap A pointer set within the fuction that points to the record
 *                   and bitmap of intrest
 *
 * @retval number of words in the bitmap
 *
 * Example Usage:
 * @code
 *      char *wah_file_name="data/10.1e4.ind.wahbm";
 *      struct wah_file wf = init_wah_file(wah_file_name);
 *      unsigned int *wah_bm;
 *      unsigned int test_record = 1;
 *      unsigned int test_bitmap = 2;
 *      unsinged int wah_size = get_wah_bitmap(wf,
 *                                             test_record,
 *                                             test_bitmap,
 *                                             &wah_bm);
 *      unsinged int *ints;
 *      unsigned int num_ints = wah_to_ints(wah_bm, wah_size, &ints);
 * @endcode
 */
unsigned int get_wah_bitmap(struct wah_file wf,
                            unsigned int wah_record,
                            unsigned int bitmap,
                            unsigned int **wah_bitmap);

/**
 * @brief Get a pointer to the bitmap of a particular WAH-encoded bitmap record.
 *
 * This is identical to get_wah_bitmap, but assumes that space for wah_bm has
 * already been allocated.
 *
 * @param wf The WAH file data structure
 * @param wah_record The record ID
 * @param bitmap The bitmap ID (0,1,2, or 3)
 * @param wah_bitmap A pointer set within the fuction that points to the record
 *                   and bitmap of intrest
 *
 * @retval number of words in the bitmap
 *
 * Example Usage:
 * @code
 *      char *wah_file_name="data/10.1e4.ind.wahbm";
 *      struct wah_file wf = init_wah_file(wah_file_name);
 *      unsigned int *wah_bm;
 *      unsigned int test_record = 1;
 *      unsigned int test_bitmap = 2;
 *      unsinged int wah_size = get_wah_bitmap_in_place(wf,
 *                                                      test_record,
 *                                                      test_bitmap,
 *                                                      &wah_bm);
 *      unsinged int *ints;
 *      unsigned int num_ints = wah_to_ints(wah_bm, wah_size, &ints);
 * @endcode
 */
unsigned int get_wah_bitmap_in_place(struct wah_file wf,
                                     unsigned int wah_record,
                                     unsigned int bitmap,
                                     unsigned int **wah_bitmap);











/**
 * @brief Get a pointer to the bitmap of a particular WAH-encoded record
 *
 * @param wf The WAH file data structure
 * @param wah_record The record ID
 * @param wah A pointer set within the fuction that points to the record
 *            of intrest
 *
 * @retval number of words in the record
 *
 * Example Usage:
 * @code
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *      struct wah_file wf = init_wah_file(wah_file_name);
 *      unsigned int *wah;
 *      unsinged int wah_size = get_wah(wf,
 *                                      test_record,
 *                                      &wah);
 *      unsinged int *ints;
 *      unsigned int num_ints = wah_to_ints(wah, wah_size, &ints);
 * @endcode
 */
unsigned int get_wah_record(struct wah_file wf,
                            unsigned int wah_record,
                            unsigned int **wah);
/**
 * @brief Get a pointer to the packed int representation of a particular plain
 * text record
 *
 * @param pf The initialized plain text file data structure
 * @param plt_record The record ID
 * @param plt A pointer set within the fuction that points to the record
 *            of intrest
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      struct plt_file pf = init_wah_file(plt_file_name);
 *      unsigned int *plt;
 *      unsinged int plt_size = get_plt_record(pf,
 *                                             test_record,
 *                                             &plt);
 * @endcode
 */
unsigned int get_plt_record(struct plt_file pf,
                            unsigned int plt_record,
                            unsigned int **plt);


/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value
 *
 * @param pf The initialized plain text encoded file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 unsigned int range_records_plt(struct plt_file pt,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int **R);
/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value using in_place functions
 *
 * @param pf The initialized plain text encoded file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */

unsigned int range_records_in_place_wahbm(struct wah_file wf,
                                          unsigned int *record_ids,
                                          unsigned int num_r,
                                          unsigned int start_test_value,
                                          unsigned int end_test_value,
                                          unsigned int **R);



/**
 * @brief Return records whose value is equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int eq_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);

/**
 * @brief Return records whose value is NOT equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int ne_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);

/**
 * @brief Return records whose value are greater than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int gt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);

/**
 * @brief Return records whose value are greater than or equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int gte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);


/**
 * @brief Return records whose value are less than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int lt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);

/**
 * @brief Return records whose value are less or equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int lte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R);


/**
 * @brief
 *
 * @param
 *
 * @returnval
 *
 * Example Usage:
 * @code
 * @endcode
 */

unsigned int gt_records_ubin(struct ubin_file uf,
                             unsigned int *record_ids,
                             unsigned int num_r,
                             unsigned int test_value,
                             unsigned int **R);


/**
 * @brief Return records whose values are >= start_test_value and < end_test_value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param start_test_value is the lower bound value to test fields against (inclusive)
 * @param end_test_value is the upper bound value to test fields against (exclusive)
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 unsigned int range_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int **R);

/**
 * @brief Return records whose values are >= start_test_value,  < end_test_value and not exclude_value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param start_test_value is the lower bound value to test fields against
 * @param end_test_value is the upper bound value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 unsigned int range_records_w_exclude_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int exclude_value,
                              unsigned int **R);

/**
 * @brief Return records whose value is equal to test_value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int eq_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);

/**
 * @brief Return records whose value is NOT equal to test_value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to exlude upon
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int ne_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);

/**
 * @brief Return records whose value are greater than the test value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int gt_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);
/**
 * @brief Return records whose value are greater than the test value using 
 * in-place functions
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int gt_records_in_place_wahbm(struct wah_file wf,
                                       unsigned int *record_ids,
                                       unsigned int num_r,
                                       unsigned int test_value,
                                       unsigned int **R);


/**
 * @brief Return records whose value are greater or equal to the test value 
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int gte_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);

/**
 * @brief Return records whose value are less than the test value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int lt_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);

/**
 * @brief Return records whose value are less than or equal to the test value 
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
unsigned int lte_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R);

/**
 * @brief 
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */

unsigned int gt_records_in_place_wahbm(struct wah_file wf,
                                       unsigned int *record_ids,
                                       unsigned int num_r,
                                       unsigned int test_value,
                                       unsigned int **R);


/**
 * @brief convert an uncompressed binary file using WAH (no bitmaps)
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 32-bit)
 * 1st record 
 * 2nd record 
 * ...
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out WAH encoded ubin file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wah(ubin_file_name, wah_file_name);
 * @endcode
 */
unsigned int convert_file_by_name_ubin_to_wah(char *ubin_in, char *wah_out);


/**
 * @brief Print a plaint text file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param pf An intitialized plaint text file
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 *
 * @retval number of records printed 
 */
unsigned int print_plt(struct plt_file pf,
                       unsigned int *record_ids,
                       unsigned int num_r);

/**
 * @brief Print plain text file (wrapper around print_plt). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param pf_file_name Plain text file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 *
 * @retval number of records printed 
 */
unsigned int print_by_name_plt(char *pf_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r);


/**
 * @brief Print an uncompressed binary file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param uf An intitialized uncompressed binary file
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *
 * @retval number of records printed 
 */
unsigned int print_ubin(struct ubin_file uf,
                        unsigned int *record_ids,
                        unsigned int num_r,
                        unsigned int format);


/**
 * @brief Print uncompressed binary file (wrapper around print_ubin). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param ubin_file_name Uncompressed binary file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *
 *
 * @retval number of records printed 
 */
unsigned int print_by_name_ubin(char *ubin_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format);
/**
 * @brief Print a WAH encoded bitmap file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param wf An initilized WAH bitmap file 
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *                              2:packed int
 *
 *
 * @retval number of records printed 
 */

unsigned int print_wahbm(struct wah_file wf,
                        unsigned int *record_ids,
                        unsigned int num_r,
                        unsigned int format);

/**
 * @brief Print a WAH encoded bitmap file (wrapper around print_wahbm). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param wahbm_file_name WAH bitmap file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *                              2:bm plain text
 *                              3:bm packed int
 *                              4:bm wah
 *
 *
 * @retval number of records printed 
 */
unsigned int print_by_name_wahbm(char *wahbm_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format);

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

unsigned int print_wah(struct wah_file wf,
                       unsigned int *record_ids,
                       unsigned int num_r,
                       unsigned int format);

/**
 * @brief Print a WAH encoded (non-bitmap) file (wrapper around print_wah). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param wahbm_file_name WAH (non-bitmap) file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *                              3:wah
 *
 *
 * @retval number of records printed 
 */
unsigned int print_by_name_wah(char *wahbm_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format);












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

void parse_cmd_line_int_csv(unsigned int *I,
                            int num_I,
                            char *cmd_line_arg);
#endif
