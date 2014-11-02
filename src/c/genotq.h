#ifndef __GENOTQ_H__
#define __GENOTQ_H__

#include <stdint.h>
#include "pq.h"
#include "bcf.h"
#include "plt.h"
#include "ubin.h"
#include "wah.h"
#include "wahbm.h"
#include "wahbm_in_place.h"
#include "wahbm_compressed_in_place.h"

struct uint_ll {
        unsigned int value;
            struct uint_ll *next;
};

unsigned int bin_char_to_int(char *bin);

int *unpack_2_bit_ints(unsigned int packed_ints);

int *unpack_1_bit_ints(unsigned int packed_ints);


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


void parse_cmd_line_int_csv(unsigned int *I,
                            int num_I,
                            char *cmd_line_arg);
const char *int_to_binary(unsigned int x);
#endif
