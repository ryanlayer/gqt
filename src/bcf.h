#ifndef __BCF_H___
#define __BCF_H___

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include "genotq.h"

struct bcf_file {
    char *file_name;
    htsFile *fp;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    uint32_t num_records;
    int32_t *gt;
};

int nlz1(unsigned x);

struct bcf_file init_bcf_file(char *file_name);

void close_bcf_file(struct bcf_file *bcf_f);

int read_unpack_next_bcf_line(struct bcf_file *bcf_f,
                              int *num_samples,
                              int *num_gts_per_sample);

uint32_t pack_sum_count_prefix_bcf_line(struct bcf_file bcf_f,
                                        uint32_t num_samples,
                                        uint32_t num_gts_per_sample,
                                        uint32_t **packed_ints,
                                        uint32_t *sum,
                                        uint32_t *prefix_len);

uint32_t md_bcf_line(struct bcf_file bcf_f,
                     char **md);

void push_bcf_gt_md(pri_queue *q,
                    struct bcf_file *bcf_f,
                    uint64_t *md_index,
                    uint32_t num_inds,
                    uint32_t num_vars,
                    char *gt_of_name,
                    char *md_of_name);


void sort_gt_md(pri_queue *q,
                uint64_t *md_index,
                uint64_t *md_lens,
                uint32_t num_inds,
                uint32_t num_vars,
                char *gt_of_name,
                char *s_gt_of_name,
                char *md_of_name,
                char *bim_of_name,
                char *vid_out);

void rotate_gt(uint32_t num_inds,
               uint32_t num_vars,
               char *s_gt_of_name,
               char *r_s_gt_of_name);

void compress_md(struct bcf_file *bcf_f,
                 char *md_of_name,
                 char *bim_out,
                 uint64_t *md_lens,
                 uint32_t num_vars);

int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out,
                                          char *vid_out,
                                          char *tmp_dir);
#endif
