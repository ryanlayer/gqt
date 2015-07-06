#ifndef __BM_H__
#define __BM_H__

/**
 * @brief Given a set of bin ranges convert an array of unsigned 32-bit ints to
 * a bitmap index
 *
 * bin_ranges defines the lower (inclusive) and upper (exclusive) range for
 * each bin and each bin corresponds to a binary array in the bitmap
 *
 * @param I array of unsigned integers to convert to a bitmap
 * @param len_I number of elements in I
 * @param bin_range_lo array defining the lower bound of each each bin
 * @param bin_range_hi array defining the upper bound of each each bin
 * @param len_bin_ranges number of bins defined in bin_range_lo and
 * bin_range_hi
 * @param less_than_bin if set to 1, the first bin will represent all values
 * less than the lowest bin range
 * @param greater_than_bin if set to 1, the last bin will represent all values
 * greater than or equal to the highest bin range
 * @param num_bins total number of bin (should equal len_bin_ranges +
 * less_than_bin + greater_than_bin)
 * @param bit_maps array containing all binary arrays, if NULL then space will
 * be allocated
 * 
 * @retval number of ints in each binary array
 *
 * Example Usage:
 * @code
 *
 *      uint32_t I[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10 };
 *
 *      uint32_t bin_range_lo_I[3] = { 3, 6, 7 };
 *      uint32_t bin_range_hi_I[3] = { 5, 7, 9 };
 *
 *      uint32_t *bm = NULL;
 *
 *      uint32_t bit_array_lens = uint32_t_to_bitmap(I,
 *                                                   10,
 *                                                   bin_range_lo_I,
 *                                                   bin_range_hi_I,
 *                                                   3,
 *                                                   1,
 *                                                   1,
 *                                                   5,
 *                                                   &bm); 
 * @endcode
 */
uint32_t uint32_t_to_bitmap(uint32_t *I,
                            uint32_t len_I,
                            uint32_t *bin_range_lo,
                            uint32_t *bin_range_hi,
                            uint32_t len_bin_ranges,
                            uint32_t less_than_bin,
                            uint32_t greater_than_bin,
                            uint32_t num_bins,
                            uint32_t **bit_map);

/**
 * @brief The double version of uint32_t_to_bitmap
 */
uint32_t double_to_bitmap(double *I,
                          uint32_t len_I,
                          double *bin_range_lo,
                          double *bin_range_hi,
                          uint32_t len_bin_ranges,
                          uint32_t less_than_bin,
                          uint32_t greater_than_bin,
                          uint32_t num_bins,
                          uint32_t **bit_map);

uint32_t int_to_bitmap(int *I,
                       uint32_t len_I,
                       int *bin_range_lo,
                       int *bin_range_hi,
                       uint32_t len_bin_ranges,
                       uint32_t less_than_bin,
                       uint32_t greater_than_bin,
                       uint32_t num_bins,
                       uint32_t **bit_map);

uint32_t float_to_bitmap(float *I,
                         uint32_t len_I,
                         float *bin_range_lo,
                         float *bin_range_hi,
                         uint32_t len_bin_ranges,
                         uint32_t less_than_bin,
                         uint32_t greater_than_bin,
                         uint32_t num_bins,
                         uint32_t **bit_map);


uint32_t int_equal_freq_binning(int *I,
                                uint32_t to_test,
                                uint32_t num_bins,
                                int **bin_range_lo,
                                int **bin_range_hi);

uint32_t float_equal_freq_binning(float *I,
                                  uint32_t to_test,
                                  uint32_t num_bins,
                                  float **bin_range_lo,
                                  float **bin_range_hi);
#endif
