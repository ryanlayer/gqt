#ifndef __VID_H__
#define __VID_H__

#include <stdint.h>
#include <stdio.h>

#include "genotq.h"

struct vid_file {
    union {
        FILE *local;
        knetFile *remote;
    } file;
    enum {
        VID_LOCAL,
        VID_REMOTE
    } type;
    char *file_name;
    uint32_t *vids;
    struct gqt_file_header *gqt_header;
};

/**
 * @brief   Open a new file and write out the GQT file header, leaving
 *          the file handle open and vids set to NULL.
 *
 * @param file_name     path to VID file
 * @param full_cmd      GQT command used to create the index, hashed for an id
 * @param num_variants  number of variants
 * @param num_samples   number of samples
 *
 * @retval pointer to a vid_file structure
 *
 * Example Usage:
 * @code
 *      char *full_cmd = "gqt convert bcf -i test.bcf";
 *      struct vid_file *v0 = new_vid_file("test.bcf.vid",
 *                                         full_cmd,
 *                                         NUM_VARS,
 *                                         NUM_INDS); 
 *      uint32_t i;
 *      for (i=0; i<NUM_VARS; ++i)
 *          write_vid(v0, (i+1));
 *
 *      destroy_vid_file(v0);
 * @endcode
 */
struct vid_file *new_vid_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples);

/**
 * @brief   Write a single uint32_t to the end of the VID file
 *
 * @param v     A vid_file that was created from new_vid_file
 * @param vid   The uint32_t that is to be appended
 *
 * Example Usage:
 * @code
 *      char *full_cmd = "gqt convert bcf -i test.bcf";
 *      struct vid_file *v0 = new_vid_file("test.bcf.vid",
 *                                         full_cmd,
 *                                         NUM_VARS,
 *                                         NUM_INDS); 
 *      uint32_t i;
 *      for (i=0; i<NUM_VARS; ++i)
 *          write_vid(v0, (i+1));
 *
 *      destroy_vid_file(v0);
 * @endcode
 */
void write_vid(struct vid_file *v, uint32_t vid);

/**
 * @brief   Open an exisiting VID file and read the header but not the data.
 *
 * @param file_name     path to VID file
 *
 * @retval pointer to a vid_file structure
 *
 * Example Usage:
 * @code
 *      struct vid_file *v1 = open_vid_file("test_vid");
 *      load_vid_data(v1);
 *
 *      for (i=0; i<NUM_VARS; ++i)
 *          printf("%d\n", v1->vids[i]);
 *
 *      destroy_vid_file(v1);
 * @endcode
 */ 
struct vid_file *open_vid_file(char *file_name);

/**
 * @brief   Close VID file and free memory
 *
 * @param vid_file  opened VID file
 *
 * Example Usage:
 * @code
 *      struct vid_file *v1 = open_vid_file("test_vid");
 *      load_vid_data(v1);
 *
 *      for (i=0; i<NUM_VARS; ++i)
 *          printf("%d\n", v1->vids[i]);
 *
 *      destroy_vid_file(v1);
 * @endcode
 */ 
void destroy_vid_file(struct vid_file *v);

/**
 * @brief   Load data from an open VID file
 *
 * @param vid_file  opened VID file
 *
 * Example Usage:
 * @code
 *      struct vid_file *v1 = open_vid_file("test_vid");
 *      load_vid_data(v1);
 *
 *      for (i=0; i<NUM_VARS; ++i)
 *          printf("%d\n", v1->vids[i]);
 *
 *      destroy_vid_file(v1);
 * @endcode
 */ 
void load_vid_data(struct vid_file *v);

#endif
