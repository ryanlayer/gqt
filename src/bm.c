#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "bm.h"
#include "genotq.h"

//{{{uint32_t uint32_t_to_bitmap(uint32_t *I,
uint32_t uint32_t_to_bitmap(uint32_t *I,
                            uint32_t len_I,
                            uint32_t *bin_range_lo,
                            uint32_t *bin_range_hi,
                            uint32_t len_bin_ranges,
                            uint32_t less_than_bin,
                            uint32_t greater_than_bin,
                            uint32_t num_bins,
                            uint32_t **bit_map)
{

    uint32_t num_ints = 1 + ((len_I - 1) / 32);

    uint32_t *target_bit_map;

    if (*bit_map == NULL) 
        *bit_map = (uint32_t *)malloc(num_ints*num_bins*sizeof(uint32_t));

    memset(*bit_map, 0, num_ints*num_bins*sizeof(uint32_t));

    uint32_t key,is_in_bin, bin, pos, i, bit_i = 0, int_i = 0;

    for (i = 0; i < len_I; ++i ) {

        key = I[i];

        pos = bsearch_uint32_t(key,
                               bin_range_lo,
                               len_bin_ranges,
                               -1,
                               len_bin_ranges);

        is_in_bin = 0;

        // first we need to check to see if the result is in either the
        // less than or great than bins
        if ( (pos == 0) && (key < bin_range_lo[pos]) ) {
            if (less_than_bin > 0) {
                is_in_bin = 1;
                bin = 0;
            } else {
                is_in_bin = 0;
            }
        } else if ( (pos == len_bin_ranges) && 
                    (key >= bin_range_hi[pos - 1]) ) {
            if (greater_than_bin > 0) {
                is_in_bin = 1;
                bin = pos + less_than_bin;
            } else {
                is_in_bin = 0;
            }
        // check if matches the lower bound
        } else if ( (pos < len_bin_ranges) && (key == bin_range_lo[pos]) ) {
            is_in_bin = 1;
            bin = pos + less_than_bin;
        // if greater than lower, check if less than upper 
        } else if ( (pos > 0) && key < bin_range_hi[pos - 1]) { 
            is_in_bin = 1;
            bin = (pos - 1) + less_than_bin;
        } 
        
        if (is_in_bin == 1) {
            // set the bit in the binary array associated with this bin
            target_bit_map = *bit_map + num_ints*bin;
            target_bit_map[int_i] |= 1 << (31 - bit_i);
        }

        bit_i += 1;

        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    return num_ints;
}
//}}}

//{{{uint32_t int_to_bitmap(int *I,
uint32_t int_to_bitmap(int *I,
                       uint32_t len_I,
                       int *bin_range_lo,
                       int *bin_range_hi,
                       uint32_t len_bin_ranges,
                       uint32_t less_than_bin,
                       uint32_t greater_than_bin,
                       uint32_t num_bins,
                       uint32_t **bit_map)
{

    uint32_t num_ints = 1 + ((len_I - 1) / 32);

    uint32_t *target_bit_map;

    if (*bit_map == NULL) 
        *bit_map = (uint32_t *)malloc(num_ints*num_bins*sizeof(uint32_t));

    memset(*bit_map, 0, num_ints*num_bins*sizeof(uint32_t));

    uint32_t key,is_in_bin, bin, pos, i, bit_i = 0, int_i = 0;

    for (i = 0; i < len_I; ++i ) {

        key = I[i];

        pos = bsearch_int(key,
                          bin_range_lo,
                          len_bin_ranges,
                          -1,
                          len_bin_ranges);

        is_in_bin = 0;

        // first we need to check to see if the result is in either the
        // less than or great than bins
        if ( (pos == 0) && (key < bin_range_lo[pos]) ) {
            if (less_than_bin > 0) {
                is_in_bin = 1;
                bin = 0;
            } else {
                is_in_bin = 0;
            }
        } else if ( (pos == len_bin_ranges) && 
                    (key >= bin_range_hi[pos - 1]) ) {
            if (greater_than_bin > 0) {
                is_in_bin = 1;
                bin = pos + less_than_bin;
            } else {
                is_in_bin = 0;
            }
        // check if matches the lower bound
        } else if ( (pos < len_bin_ranges) && (key == bin_range_lo[pos]) ) {
            is_in_bin = 1;
            bin = pos + less_than_bin;
        // if greater than lower, check if less than upper 
        } else if ( (pos > 0) && key < bin_range_hi[pos - 1]) { 
            is_in_bin = 1;
            bin = (pos - 1) + less_than_bin;
        } 
        
        if (is_in_bin == 1) {
            // set the bit in the binary array associated with this bin
            target_bit_map = *bit_map + num_ints*bin;
            target_bit_map[int_i] |= 1 << (31 - bit_i);
        }

        bit_i += 1;

        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    return num_ints;
}
//}}}

//{{{uint32_t double_to_bitmap(double *I,
uint32_t double_to_bitmap(double *I,
                          uint32_t len_I,
                          double *bin_range_lo,
                          double *bin_range_hi,
                          uint32_t len_bin_ranges,
                          uint32_t less_than_bin,
                          uint32_t greater_than_bin,
                          uint32_t num_bins,
                          uint32_t **bit_map)
{
    uint32_t num_ints = 1 + ((len_I - 1) / 32);

    uint32_t *target_bit_map;

    if (*bit_map == NULL) 
        *bit_map = (uint32_t *)malloc(num_ints*num_bins*sizeof(uint32_t));

    memset(*bit_map, 0, num_ints*num_bins*sizeof(uint32_t));

    uint32_t is_in_bin, bin, pos, i, bit_i = 0, int_i = 0;
    double key;

    for (i = 0; i < len_I; ++i ) {

        key = I[i];

        pos = bsearch_double(key,
                             bin_range_lo,
                             len_bin_ranges,
                             -1,
                             len_bin_ranges);

        is_in_bin = 0;

        // first we need to check to see if the result is in either the
        // less than or great than bins
        if ( (pos == 0) && (key < bin_range_lo[pos]) ) {
            if (less_than_bin > 0) {
                is_in_bin = 1;
                bin = 0;
            } else {
                is_in_bin = 0;
            }
        } else if ( (pos == len_bin_ranges) && 
                    (key >= bin_range_hi[pos - 1]) ) {
            if (greater_than_bin > 0) {
                is_in_bin = 1;
                bin = pos + less_than_bin;
            } else {
                is_in_bin = 0;
            }
        // check if matches the lower bound
        } else if ( (pos < len_bin_ranges) && (key == bin_range_lo[pos]) ) {
            is_in_bin = 1;
            bin = pos + less_than_bin;
        // if greater than lower, check if less than upper 
        } else if ( (pos > 0) && key < bin_range_hi[pos - 1]) { 
            is_in_bin = 1;
            bin = (pos - 1) + less_than_bin;
        } 
        
        if (is_in_bin == 1) {
            // set the bit in the binary array associated with this bin
            target_bit_map = *bit_map + num_ints*bin;
            target_bit_map[int_i] |= 1 << (31 - bit_i);
        }

        bit_i += 1;

        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    return num_ints;
}
//}}}

//{{{ uint32_t float_to_bitmap(float *I,
uint32_t float_to_bitmap(float *I,
                         uint32_t len_I,
                         float *bin_range_lo,
                         float *bin_range_hi,
                         uint32_t len_bin_ranges,
                         uint32_t less_than_bin,
                         uint32_t greater_than_bin,
                         uint32_t num_bins,
                         uint32_t **bit_map)
{
    uint32_t num_ints = 1 + ((len_I - 1) / 32);

    uint32_t *target_bit_map;

    if (*bit_map == NULL) 
        *bit_map = (uint32_t *)malloc(num_ints*num_bins*sizeof(uint32_t));

    memset(*bit_map, 0, num_ints*num_bins*sizeof(uint32_t));

    uint32_t is_in_bin, bin, pos, i, bit_i = 0, int_i = 0;
    float key;

    for (i = 0; i < len_I; ++i ) {

        key = I[i];

        pos = bsearch_float(key,
                            bin_range_lo,
                            len_bin_ranges,
                            -1,
                            len_bin_ranges);

        is_in_bin = 0;

        // first we need to check to see if the result is in either the
        // less than or great than bins
        if ( (pos == 0) && (key < bin_range_lo[pos]) ) {
            if (less_than_bin > 0) {
                is_in_bin = 1;
                bin = 0;
            } else {
                is_in_bin = 0;
            }
        } else if ( (pos == len_bin_ranges) && 
                    (key >= bin_range_hi[pos - 1]) ) {
            if (greater_than_bin > 0) {
                is_in_bin = 1;
                bin = pos + less_than_bin;
            } else {
                is_in_bin = 0;
            }
        // check if matches the lower bound
        } else if ( (pos < len_bin_ranges) && (key == bin_range_lo[pos]) ) {
            is_in_bin = 1;
            bin = pos + less_than_bin;
        // if greater than lower, check if less than upper 
        } else if ( (pos > 0) && key < bin_range_hi[pos - 1]) { 
            is_in_bin = 1;
            bin = (pos - 1) + less_than_bin;
        } 
        
        if (is_in_bin == 1) {
            // set the bit in the binary array associated with this bin
            target_bit_map = *bit_map + num_ints*bin;
            target_bit_map[int_i] |= 1 << (31 - bit_i);
        }

        bit_i += 1;

        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    return num_ints;
}
//}}}

//{{{ uint32_t int_equal_freq_binning(int *I,
uint32_t int_equal_freq_binning(int *I,
                                uint32_t to_test,
                                uint32_t num_bins,
                                int **bin_range_lo,
                                int **bin_range_hi)
{

    int *sorted_I = (int *)malloc(to_test * sizeof(int));
    memcpy(sorted_I, I, to_test * sizeof(int));

    qsort(sorted_I, to_test, sizeof(int), int_compare);

    uint32_t i, num_uniq = 1;

    // create a list of uniq values on top of the sorted_I
    for (i=1; i < to_test; ++i) {
        if (sorted_I[i-1] != sorted_I[i]) {
            sorted_I[num_uniq] = sorted_I[i];
            num_uniq += 1;
        }
    }

    if (num_uniq < num_bins)
        num_bins = num_uniq;

    *bin_range_lo = (int *)malloc(num_bins * sizeof(int));
    *bin_range_hi = (int *)malloc(num_bins * sizeof(int));

    float rough_bin_freq = ((float)num_uniq) / ((float)num_bins);

    (*bin_range_lo)[0] = sorted_I[0];
    (*bin_range_hi)[num_bins-1] = sorted_I[num_uniq-1] + 1;

    for (i=1; i < num_bins; ++i) {
        (*bin_range_hi)[i-1] = sorted_I[(int)(i*rough_bin_freq)];
        (*bin_range_lo)[i] = sorted_I[(int)(i*rough_bin_freq)];
    }

    return num_bins;
}
//}}}

//{{{ uint32_t float_equal_freq_binning(float *I,
uint32_t float_equal_freq_binning(float *I,
                                  uint32_t to_test,
                                  uint32_t num_bins,
                                  float **bin_range_lo,
                                  float **bin_range_hi)
{

    float *sorted_I = (float *)malloc(to_test * sizeof(float));
    memcpy(sorted_I, I, to_test * sizeof(float));

    uint32_t i, num_uniq = 1;

    qsort(sorted_I, to_test, sizeof(float), float_compare);

    // create a list of uniq values on top of the sorted_I
    for (i=1; i < to_test; ++i) {
        if (sorted_I[i-1] != sorted_I[i]) {
            sorted_I[num_uniq] = sorted_I[i];
            num_uniq += 1;
        }
    }

    if (num_uniq < num_bins)
        num_bins = num_uniq;

    *bin_range_lo = (float *)malloc(num_bins * sizeof(float));
    *bin_range_hi = (float *)malloc(num_bins * sizeof(float));

    float rough_bin_freq = ((float)num_uniq) / ((float)num_bins);

    (*bin_range_lo)[0] = sorted_I[0];
    (*bin_range_hi)[num_bins-1] = sorted_I[num_uniq-1] + 1;

    for (i=1; i < num_bins; ++i) {
        (*bin_range_hi)[i-1] = sorted_I[(int)(i*rough_bin_freq)];
        (*bin_range_lo)[i] = sorted_I[(int)(i*rough_bin_freq)];
    }

    return num_bins;
}
//}}}
