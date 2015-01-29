#ifndef __WAHBM_COMPRESSED_IN_PLACE_H__
#define __WAHBM_COMPRESSED_IN_PLACE_H__

#include <stdint.h>
#include "genotq.h"

/**
 * @brief   compressed in-place OR two WAH runs
 *
 * The assumption here is that r_wah is pre-allocated to an array that can hold
 * the maximum possible size. r_wah can contain a mix of fills and literals,
 * but any fill must be followed by enough empty words to allow the fill to be
 * completely inflated if future operations require it.
 * wah can be any combination of litterals and fills.  The result is stored
 * back into r_wah.
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
 *    uint32_t X[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *        };
 *
 *    uint32_t Y[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *        };
 *    uint32_t R[6] = {0x80000006,0,0,0,0,0};
 *
 *    uint32_t r = wah_compressed_in_place_or(R, 6, w_Y, wah_size_Y);
 *    r = wah_compressed_in_place_or(R, 6, w_X, wah_size_X);
 * @endcode
 */
uint32_t  wah_compressed_in_place_or(uint32_t *r_wah,
                                         uint32_t r_wah_size,
                                         uint32_t *wah,
                                         uint32_t wah_size);
/**
 * @brief   compressed in-place AND two WAH runs
 *
 * The assumption here is that r_wah is pre-allocated to an array that can hold
 * the maximum possible size. r_wah can contain a mix of fills and literals,
 * but any fill must be followed by enough empty words to allow the fill to be
 * completely inflated if future operations require it.
 * wah can be any combination of litterals and fills.  The result is stored
 * back into r_wah.
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
 *    uint32_t X[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("01000000000101010100000000000000"),
 *          bin_char_to_int("01000000000000000001010101000000")
 *        };
 *
 *    uint32_t Y[5] =
 *        { bin_char_to_int("01000000000000000000000000000001"),
 *          bin_char_to_int("11111111111111111111111111111111"),
 *          bin_char_to_int("11111111111111111111111111111000"),
 *          bin_char_to_int("00000000000000000000000000000000"),
 *          bin_char_to_int("00000000000000000000000000001011")
 *        };
 *    uint32_t R[6] = {0xc0000006,0,0,0,0,0};
 *
 *    uint32_t r = wah_compressed_in_place_and(R, 6, w_Y, wah_size_Y);
 *    r = wah_compressed_in_place_and(R, 6, w_X, wah_size_X);
 * @endcode
 */
uint32_t  wah_compressed_in_place_and(uint32_t *r_wah,
                                          uint32_t r_wah_size,
                                          uint32_t *wah,
                                          uint32_t wah_size);

/**
 * @brief   AND two compressed in-place WAH runs
 *
 * The assumption here is that r_wah and whar are pre-allocated to arrays
 * that can hold the maximum possible size. Both can contain a mix of fills
 * and literals, but any fill must be followed by enough empty words to allow
 * the fill to be completely inflated if future operations require it.  
 *
 * @param r_wah a compressed in-place wah-encoded bitmap of literals and fills
 * @param r_wah_size the number of words in r_wah
 * @param wah a compressed in-place wah-encoded bitmap of literals and fills
 * @param wah_size the numberof words in wah
 *
 * @retval  The number of elements in r_wah
 *
 * @ingroup WAH
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t  wah_compressed_in_place_and_compressed_in_place(
                                          uint32_t *r_wah,
                                          uint32_t r_wah_size,
                                          uint32_t *wah,
                                          uint32_t wah_size);


/**
 * @brief Convert compressed in place WAH encoding to uncompressed binary ints.
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

uint32_t compressed_in_place_wah_to_ints(uint32_t *W,
                                             uint32_t W_len,
                                             uint32_t **O);
/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value using compressed in place functions
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
uint32_t range_records_compressed_in_place_wahbm(
            struct wah_file wf,
            uint32_t *record_ids,
            uint32_t num_r,
            uint32_t start_test_value,
            uint32_t end_test_value,
            uint32_t **R);


/**
 * @brief For each field, count that number of records that meet a certian
 * critera using compressed in-place functions
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param start_test_value is the lower bound value to test fields against
 * (inclusive)
 * @param end_test_value is the upper bound value to test fields against
 * (exclusive)
 * @param R interger counts
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t count_range_records_compressed_in_place_wahbm(
            struct wah_file wf,
            uint32_t *record_ids,
            uint32_t num_r,
            uint32_t start_test_value,
            uint32_t end_test_value,
            uint32_t **R);


/**
 * @brief Take the numbers encoded by a single compressed in-place wahbm, and
 * add them to the R and store back into R
 *
 * @param R array of ints that are added to the WAH-encoded bit map store 
 * @param r_size size of R
 * @param wah A compressed in-place WAH-encoded bitmap
 * @param wah_size number of ints in wah
 *
 * @return the size of R
 *
 */
uint32_t add_compressed_in_place_wahbm(uint32_t *R,
                                           uint32_t r_size,
                                           uint32_t *wah,
                                           uint32_t wah_size);

/**
 * @brief For each field, count the number of records that meet the criterea
 * using compressed in-place funcitons
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R integer counts
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gt_count_records_compressed_in_place_wahbm(
                    struct wah_file wf,
                    uint32_t *record_ids,
                    uint32_t num_r,
                    uint32_t test_value,
                    uint32_t **R);

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
uint32_t gt_records_compressed_in_place_wahbm(struct wah_file wf,
                                       uint32_t *record_ids,
                                       uint32_t num_r,
                                       uint32_t test_value,
                                       uint32_t **R);

uint32_t avx_sum_range_records_in_place_wahbm(
            struct wah_file wf,
            uint32_t *record_ids,
            uint32_t num_r,
            uint32_t start_test_value,
            uint32_t end_test_value,
            uint32_t **R);

#endif
