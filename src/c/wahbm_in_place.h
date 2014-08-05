#ifndef __WAHBM_IN_PLACE_H__
#define __WAHBM_IN_PLACE_H__

#include <stdint.h>
#include "genotq.h"

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
unsigned int count_range_records_in_place_wahbm(struct wah_file wf,
                                                unsigned int *record_ids,
                                                unsigned int num_r,
                                                unsigned int start_test_value,
                                                unsigned int end_test_value,
                                                unsigned int **R);

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
unsigned int gt_count_records_in_place_wahbm(struct wah_file wf,
                                             unsigned int *record_ids,
                                             unsigned int num_r,
                                             unsigned int test_value,
                                             unsigned int **R) ;

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
unsigned int sum_range_records_in_place_wahbm(struct wah_file wf,
                                              unsigned int *record_ids,
                                              unsigned int num_r,
                                              unsigned int start_test_value,
                                              unsigned int end_test_value,
                                              unsigned int **R);


unsigned int gt_sum_records_in_place_wahbm(struct wah_file wf,
                                             unsigned int *record_ids,
                                             unsigned int num_r,
                                             unsigned int test_value,
                                             unsigned int **R);

#endif
