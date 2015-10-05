#ifndef __BCF_H___
#define __BCF_H___

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include "pq.h"
#include "off.h"

struct bcf_file {
    char *file_name;
    uint32_t is_bcf;
    union {
        htsFile *bcf;
        BGZF *vcf;
    } fp;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    uint32_t num_records;
    uint64_t offset;
    kstring_t str;
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

void push_bcf_gt_offset(pri_queue *q,
                       struct bcf_file *bcf_f,
                       uint32_t num_inds,
                       uint32_t num_vars,
                       char *gt_of_name,
                       char *offset_of_name,
                       char *full_cmd);

void sort_gt(pri_queue *q,
             uint32_t num_inds,
             uint32_t num_vars,
             char *gt_of_name,
             char *gt_s_of_name,
             char *vid_out,
             char *full_cmd);

void rotate_gt(uint32_t num_inds,
               uint32_t num_vars,
               char *s_gt_of_name,
               char *r_s_gt_of_name);

void compress_md(struct bcf_file *bcf_f,
                 char *md_of_name,
                 char *bim_out,
                 uint64_t *md_lens,
                 uint32_t num_vars,
                 uint32_t num_inds,
                 char *full_cmd);

int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out,
                                          char *vid_out,
                                          char *tmp_dir,
                                          char *full_cmd);

int convert_file_by_name_bcf_to_wahbm_offset(char *in,
                                             uint32_t num_fields,
                                             uint32_t num_records,
                                             char *wah_out,
                                             char *offset_out,
                                             char *vid_out,
                                             char *tmp_dir,
                                             char *full_cmd);

int get_bcf_line(struct bcf_file *bcf_f);

int goto_bcf_line(struct bcf_file *bcf_f,
                  struct off_file *off_f,
                  uint32_t line_no);
#endif
