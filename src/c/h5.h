#ifndef __H5_H__
#define __H5_H__

#include <hdf5.h>
#include "genotq.h"
#include "pq.h"

struct hdf5_file
{
    hid_t file_id,
          plist_id,
          gt_dataspace_id,
          md_dataspace_id,
          r_gt_dataspace_id;
    uint32_t num_gt_ints,
             num_r_gt_ints,
             num_vars,
             num_inds;
};

struct hdf5_file init_hdf5_file(char *hdf5_file_name,
                                uint32_t num_vars,
                                uint32_t num_inds);

int write_hdf5_gt(struct hdf5_file hdf5_f,
                  uint32_t id,
                  uint32_t *gt,
                  char *md);


int read_hdf5_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *gt);

int read_hdf5_md(struct hdf5_file hdf5_f, uint32_t id, char **md);

int init_r_gt(struct hdf5_file hdf5_f);

int set_r_gt(struct hdf5_file hdf5_f,
             uint32_t r_gt_i,
             uint32_t i,
             uint32_t v);

int read_hdf5_r_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *r_gt);

int close_hdf5_file(struct hdf5_file hdf5_f);

void sort_rotate_gt_md(pri_queue *q,
                       struct hdf5_file *hdf5_f,
                       char *md_out_file_name);
#endif
