#ifndef __BCF_H___
#define __BCF_H___

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include "genotq.h"

struct bcf_file {
    htsFile *fp;
    bcf_hdr_t *hdr;
    bcf1_t *line;
    unsigned int num_records;
    int32_t *gt;
};

int nlz1(unsigned x);

struct bcf_file init_bcf_file(char *file_name);

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
                    struct hdf5_file *hdf5_f);

int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out);
#endif
