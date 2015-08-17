#ifndef __WAHBM_IN_PLACE_H__
#define __WAHBM_IN_PLACE_H__

#include <stdint.h>
#include "wah.h"
#include "wahbm.h"

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
 *    uint32_t R[6] = {0,0,0,0,0,0};
 *
 *    uint32_t r = wah_in_place_or(R, 6, w_Y, wah_size_Y);
 *    r = wah_in_place_or(R, 6, w_X, wah_size_X);
 *    uint32_t *ints_ip;
 *    uint32_t ints_ip_len = wah_to_ints(R,r,&ints_ip);
 * @endcode
 */
uint32_t  wah_in_place_or(uint32_t *r_wah,
                              uint32_t r_wah_size,
                              uint32_t *wah,
                              uint32_t wah_size);

uint32_t wah_in_place_xor(uint32_t *r_wah,
                          uint32_t r_wah_size,
                          uint32_t *wah,
                          uint32_t wah_size);


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
 *    uint32_t R[6] = {0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff,
 *                         0x7fffffff};
 *
 *    uint32_t r = wah_in_place_and(R, 6, w_Y, wah_size_Y);
 *    r = wah_in_place_and(R, 6, w_X, wah_size_X);
 *    uint32_t *ints_ip;
 *    uint32_t ints_ip_len = wah_to_ints(R,r,&ints_ip);
 * @endcode
 */
uint32_t  wah_in_place_and(uint32_t *r_wah,
                           uint32_t r_wah_size,
                           uint32_t *wah,
                           uint32_t wah_size);

#if 0
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
 *      uint32_t *wah_bm;
 *      uint32_t test_record = 1;
 *      uint32_t test_bitmap = 2;
 *      unsinged int wah_size = get_wah_bitmap_in_place(wf,
 *                                                      test_record,
 *                                                      test_bitmap,
 *                                                      &wah_bm);
 *      unsinged int *ints;
 *      uint32_t num_ints = wah_to_ints(wah_bm, wah_size, &ints);
 * @endcode
 */
uint32_t get_wah_bitmap_in_place(struct wah_file wf,
                                 uint32_t wah_record,
                                 uint32_t bitmap,
                                 uint32_t **wah_bitmap);
#endif

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
uint32_t range_records_in_place_wahbm(struct wahbm_file *wf,
                                      uint32_t *record_ids,
                                      uint32_t num_r,
                                      uint32_t start_test_value,
                                      uint32_t end_test_value,
                                      uint32_t **R);
/**
 * @brief For each field, count that number of records that meet a certian
 * critera using in-place functions
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
uint32_t count_range_records_in_place_wahbm(struct wahbm_file *wf,
                                            uint32_t *record_ids,
                                            uint32_t num_r,
                                            uint32_t start_test_value,
                                            uint32_t end_test_value,
                                            uint32_t **R);

/**
 * @brief For each field, count the number of records that meet the criterea
 * using in-place funcitons
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
uint32_t gt_count_records_in_place_wahbm(struct wahbm_file *wf,
                                         uint32_t *record_ids,
                                         uint32_t num_r,
                                         uint32_t test_value,
                                         uint32_t **R) ;

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
uint32_t gt_records_in_place_wahbm(struct wahbm_file *wf,
                                   uint32_t *record_ids,
                                   uint32_t num_r,
                                   uint32_t test_value,
                                   uint32_t **R);

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

uint32_t gt_records_in_place_wahbm(struct wahbm_file *wf,
                                   uint32_t *record_ids,
                                   uint32_t num_r,
                                   uint32_t test_value,
                                   uint32_t **R);

/**
 * @brief For each field, sum the number of that meet a certian
 * critera using in-place functions
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
uint32_t sum_range_records_in_place_wahbm(struct wahbm_file *wf,
                                          uint32_t *record_ids,
                                          uint32_t num_r,
                                          uint32_t start_test_value,
                                          uint32_t end_test_value,
                                          uint32_t **R);


uint32_t gt_sum_records_in_place_wahbm(struct wahbm_file *wf,
                                       uint32_t *record_ids,
                                       uint32_t num_r,
                                       uint32_t test_value,
                                       uint32_t **R);
#if 0
uint32_t get_wah_bitmaps_in_place(struct wah_file wf,
                                  uint32_t wah_record,
                                  uint32_t **wah_bitmap,
                                  uint32_t *wah_sizes);
#endif

#ifdef __AVX2__
uint32_t avx_gt_count_records_in_place_wahbm(struct wahbm_file *wf,
                                             uint32_t *record_ids,
                                             uint32_t num_r,
                                             uint32_t test_value,
                                             uint32_t **R);
#endif

#ifdef __AVX2__
uint32_t avx_count_range_records_in_place_wahbm(struct wahbm_file *wf,
                                                uint32_t *record_ids,
                                                uint32_t num_r,
                                                uint32_t start_test_value,
                                                uint32_t end_test_value,
                                                uint32_t **R);
#endif

#ifdef __AVX2__
uint32_t avx_sum_range_records_in_place_wahbm(struct wahbm_file *wf,
                                              uint32_t *record_ids,
                                              uint32_t num_r,
                                              uint32_t start_test_value,
                                              uint32_t end_test_value,
                                              uint32_t **R);
#endif


#endif
