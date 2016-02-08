#ifndef __GENOTQ_H__
#define __GENOTQ_H__

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <err.h>
#include <sysexits.h>
#include <sys/errno.h>
/*
#include "parse_q.h"
#include "pq.h"
#include "ped.h"
#include "bcf.h"
#include "plt.h"
#include "ubin.h"
#include "wah.h"
#include "wahbm.h"
#include "wahbm_in_place.h"
#include "wahbm_compressed_in_place.h"
*/

#define CHUNK 16384

#define PROGRAM_NAME  "gqt"
#define MAJOR_VERSION "1"
#define MINOR_VERSION "1"
#define REVISION_VERSION "3"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION
#define MORE_SIZE 20

struct gqt_file_header {
    char marker[3]; // "GQT"
    char type; // g gqt, v vid, b bim, o off
    uint32_t major, minor, revision, build;
    uint32_t magic; //0x11223344
    unsigned long id_hash;// used to varify the files were created together.
    uint32_t num_variants, num_samples;
    uint32_t more[20];
};

struct gqt_file_header *new_gqt_file_header(char type,
                                            char *full_cmd,
                                            uint32_t num_variants,
                                            uint32_t num_samples);

struct gqt_file_header *read_gqt_file_header(char *file_name, FILE *f);

struct uint_ll {
    uint32_t value;
    struct uint_ll *next;
};

uint32_t bin_char_to_int(char *bin);

int *unpack_2_bit_ints(uint32_t packed_ints);

int *unpack_1_bit_ints(uint32_t packed_ints);

int popcount(uint32_t x);

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
uint32_t ints_to_rle(uint32_t *I, int I_len, uint32_t **O);


void parse_cmd_line_int_csv(uint32_t *I,
                            int num_I,
                            char *cmd_line_arg);
const char *int_to_binary(uint32_t x);

void check_file_read(char *file_name, FILE *fp, size_t exp, size_t obs);

int check_field_name(char *field_name);

int is_int(char *s, int *v);

unsigned long hash_cmd(char *full_cmd);
#endif
