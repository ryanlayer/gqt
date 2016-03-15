#ifndef __BIM_H__
#define __BIM_H__

#include <stdint.h>
#include <stdio.h>
#include <htslib/knetfile.h>

#include "genotq.h"

 /* 
 * The file is :
 * uncompressed size     ( sizeof(uint64_t))
 * compressed size       ( sizeof(uint64_t))
 * header size           ( sizeof(uint64_t))
 * md line lengths       ( bcf_f->num_records*sizeof(uint64_t))
 * compressed data 
 */
struct bim_file_header {
    uint64_t u_size, c_size, h_size; 
    uint64_t *md_line_lens;
};

struct bim_file {
    //FILE *file;
    union {
        FILE *local;
        knetFile *remote;
    } file;
    enum {
        BIM_LOCAL,
        BIM_REMOTE
    } type;
    char *file_name;
    uint64_t data_start;
    struct gqt_file_header *gqt_header;
    struct bim_file_header *bim_header;
};

/**
 * @brief Read the BIM specific header, which occurs after the gqt header
 *
 * This function assumes there is a GQT heeader first, and will seek to correct
 * positoin in the file
 *
 * @param bim_file An open bim_file structure
 *
 * @retval A bim_file_header with all fields, including the md lens 
 *
 * Example Usage:
 * @code
 *      struct bim_file *b = (struct bim_file *) 
 *              malloc(sizeof(struct bim_file));
 *      b->file = fopen(file_name,"rb");
 *      b->bim_header = read_bim_file_header(b);
 * @endcode
 */
struct bim_file_header *read_bim_file_header(struct bim_file *b);

/**
 * @brief Create a new BIM specific header structure
 *
 * @param u_size size of the uncompressed data
 * @param c_size size of the compressed data
 * @param h_size size of the header
 * @param md_line_lens array of line lengths
 *
 * @retval A bim_file_header with all fields, including the md lens 
 *
 * Example Usage:
 * @code
 *      uint64_t u_size = 1;
 *      uint64_t c_size = 2;
 *      uint64_t h_size = 3;
 *      uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
 *      md_line_lens[0] = 1;
 *      md_line_lens[1] = 2;
 *      md_line_lens[2] = 3;
 *      md_line_lens[3] = 4;
 *      b->bim_header = new_bim_file_header(u_size,
 *                                          c_size,
 *                                          h_size,
 *                                          md_line_lens);
 * @endcode
 */
struct bim_file_header *new_bim_file_header(uint64_t u_size,
                                            uint64_t c_size,
                                            uint64_t h_size,
                                            uint64_t *md_line_lens);

/**
 * @brief Create a new bim file and write the GQT and BIM header to file
 *
 * @param 
 * @param file_name bim file path
 * @param full_cmd command used to generate indicies
 * @param num_variants number of variants
 * @param num_samples number of samples
 * @param u_size size of the uncompressed data
 * @param c_size size of the compressed data
 * @param h_size size of the header
 * @param md_line_lens array of line lengths

 * @retval A bim_file with all fields, including the md lens 
 *
 * Example Usage:
 * @code
 *      char *file_name = "test_bim";
 *      char *full_cmd = "gqt convert bcf -i bcf";
 *      uint32_t num_variants = 4;
 *      uint32_t num_samples = 4;
 *      uint64_t u_size = 1;
 *      uint64_t c_size = 2;
 *      uint64_t h_size = 3;
 *      uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
 *      md_line_lens[0] = 1;
 *      md_line_lens[1] = 2;
 *      md_line_lens[2] = 3;
 *      md_line_lens[3] = 4;
 *
 *      struct bim_file *b = new_bim_file(file_name,
 *                                        full_cmd,
 *                                        num_variants,
 *                                        num_samples,
 *                                        u_size,
 *                                        c_size,
 *                                        h_size,
 *                                        md_line_lens);
 *
 *      destroy_bim_file(b);
 * @endcode
 */
struct bim_file *new_bim_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples,
                              uint64_t u_size,
                              uint64_t c_size,
                              uint64_t h_size,
                              uint64_t *md_line_lens);

/**
 * @brief Open a BIM file and read in GQT and BIM headers
 *
 * @param file_name path to BIM file
 *
 * @retval A bim_file with all fields, including the md lens 
 *
 * Example Usage:
 * @code
 *      struct bim_file *b = open_bim_file(file_name);
 *      destroy_bim_file(b);
 * @endcode
 */
struct bim_file *open_bim_file(char *file_name);

/**
 * @brief Update three scalors in the BIM header in struct and file
 *
 * @param u_size size of the uncompressed data
 * @param c_size size of the compressed data
 * @param h_size size of the header
 * @param b an exisiting bim_file_header
 *
 * Example Usage:
 * @code
 *      struct bim_file *b = open_bim_file(file_name);
 *      update_bim_file_header(u_size, c_size, h_size, b);
 *      destroy_bim_file(b);
 * @endcode
 */
void update_bim_file_header(uint64_t u_size,
                            uint64_t c_size,
                            uint64_t h_size,
                            struct bim_file *b);

/**
 * @brief Close BIM file and free memory in GQT and BIM headers
 *
 * @param b bim_file 
 *
 * Example Usage:
 * @code
 *      struct bim_file *b = open_bim_file(file_name);
 *      destroy_bim_file(b);
 * @endcode
 */
void destroy_bim_file(struct bim_file *b);

/**
 * @brief Jump to the start of the data in the BIM file
 *
 * @param b bim_file 
 *
 * Example Usage:
 * @code
 *      struct bim_file *b = open_bim_file(file_name);
 *      seek_bim_to_data(b)
 *      destroy_bim_file(b);
 * @endcode
 */
void seek_bim_to_data(struct bim_file *b);
#endif
