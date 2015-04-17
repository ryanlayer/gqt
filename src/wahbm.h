#ifndef __WAHBM_H__
#define __WAHBM_H__

#include <stdint.h>
#include <immintrin.h>
#include "plt.h"
#include "ubin.h"
#include "wah.h"
#include "pthread_pool.h"

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

uint32_t print_wahbm(struct wah_file wf,
                     uint32_t *record_ids,
                     uint32_t num_r,
                     uint32_t format);
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
uint32_t print_by_name_wahbm(char *wahbm_file_name,
                             uint32_t *record_ids,
                             uint32_t num_r,
                             uint32_t format);

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
 *      uint32_t *wah_bm;
 *      uint32_t test_record = 1;
 *      uint32_t test_bitmap = 2;
 *      unsinged int wah_size = get_wah_bitmap(wf,
 *                                             test_record,
 *                                             test_bitmap,
 *                                             &wah_bm);
 *      unsinged int *ints;
 *      uint32_t num_ints = wah_to_ints(wah_bm, wah_size, &ints);
 * @endcode
 */
uint32_t get_wah_bitmap(struct wah_file wf,
                        uint32_t wah_record,
                        uint32_t bitmap,
                        uint32_t **wah_bitmap);

/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value
 *
 * @param wf The initialized WAH-encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param start_test_value is the lower bound value to test fields against
 * (inclusive)
 * @param end_test_value is the upper bound value to test fields against
 * (exclusive)
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 uint32_t range_records_wahbm(struct wah_file wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t start_test_value,
                              uint32_t end_test_value,
                              uint32_t **R);

/**
 * @brief For each field, count that number of records that meet a certian
 * critera 
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

uint32_t count_range_records_wahbm(struct wah_file wf,
                                   uint32_t *record_ids,
                                   uint32_t num_r,
                                   uint32_t start_test_value,
                                   uint32_t end_test_value,
                                   uint32_t **R);

/**
 * @brief Take the numbers encoded by a single wahbm, and add them to the R
 * and store back into R
 *
 * @param R array of ints that are added to the WAH-encoded bit map store 
 * @param r_size size of R
 * @param wah A WAH-encoded bitmap
 * @param wah_size number of ints in wah
 *
 * @return the size of R
 *
 */
uint32_t add_wahbm(uint32_t *R,
                   uint32_t r_size,
                   uint32_t *wah,
                   uint32_t wah_size);

#ifdef __AVX2__
uint32_t avx_add_wahbm(uint32_t *R,
                       uint32_t r_size,
                       uint32_t *wah,
                       uint32_t wah_size);
#endif

#ifdef __AVX2__
void avx_add(uint32_t bits,
             __m256i *s_1,
             __m256i *s_2,
             __m256i *s_3,
             __m256i *s_4,
             __m256i *m,
             __m256i *R_avx,
             uint32_t avx_i);
#endif


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
uint32_t eq_records_wahbm(struct wah_file wf,
                          uint32_t *record_ids,
                          uint32_t num_r,
                          uint32_t test_value,
                          uint32_t **R);

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
uint32_t ne_records_wahbm(struct wah_file wf,
                          uint32_t *record_ids,
                          uint32_t num_r,
                          uint32_t test_value,
                          uint32_t **R);

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
uint32_t gt_records_wahbm(struct wah_file wf,
                          uint32_t *record_ids,
                          uint32_t num_r,
                          uint32_t test_value,
                          uint32_t **R);

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
uint32_t gte_records_wahbm(struct wah_file wf,
                           uint32_t *record_ids,
                           uint32_t num_r,
                           uint32_t test_value,
                           uint32_t **R);

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
uint32_t lt_records_wahbm(struct wah_file wf,
                          uint32_t *record_ids,
                          uint32_t num_r,
                          uint32_t test_value,
                          uint32_t **R);

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
uint32_t lte_records_wahbm(struct wah_file wf,
                           uint32_t *record_ids,
                           uint32_t num_r,
                           uint32_t test_value,
                           uint32_t **R);

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
 uint32_t range_records_w_exclude_wahbm(struct wah_file wf,
                                        uint32_t *record_ids,
                                        uint32_t num_r,
                                        uint32_t start_test_value,
                                        uint32_t end_test_value,
                                        uint32_t exclude_value,
                                        uint32_t **R);

/**
 * @brief For each field, count the number of records that meet the criterea
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
uint32_t gt_count_records_wahbm(struct wah_file wf,
                                uint32_t *record_ids,
                                uint32_t num_r,
                                uint32_t test_value,
                                uint32_t **R);

uint32_t add_n_wahbm(uint32_t *R,
                     uint32_t n,
                     uint32_t r_size,
                     uint32_t *wah,
                     uint32_t wah_size);

uint32_t p_pool_add_n_wahbm(uint32_t *R,
                            uint32_t n,
                            uint32_t r_size,
                            uint32_t *wah,
                            uint32_t wah_size,
                            struct pool *t_pool);

void *t_add_n_wahbm(void *arg);
void *t_add_n_wahbm_2(void *arg);

struct t_add_n_wahbm_args {
    uint32_t bits, field_i, *R, r_size, n, num_words;
};

uint32_t avx_add_n_wahbm(uint32_t *R,
                         uint32_t n,
                         uint32_t r_size,
                         uint32_t *wah,
                         uint32_t wah_size);
#ifdef __AVX2__
void avx_add_n(uint32_t bits,
               __m256i *s_1,
               __m256i *s_2,
               __m256i *s_3,
               __m256i *s_4,
               __m256i *m,
               __m256i *N,
               __m256i *R_avx,
               uint32_t field_i);
#endif

uint32_t wahbm_pca_by_name(char *in, char *out);

uint32_t wahbm_speed_check(char *in);

uint32_t wahbm_top_n_matches_by_name(char *in, uint32_t n);

uint32_t wahbm_hamm_dist_by_name(char *in, char *out);
////////////END//////////////////////////////////////////////////////
#endif
