#include <stdio.h>
#include <stdlib.h>
#include "parse_q.h"
#include "genotq.h"
#include "unity.h"
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

char *BCF_FILE = "../data/diff_gts.bcf";
uint32_t NUM_INDS = 10;
uint32_t NUM_VARS = 2;

void setUp(void) { }
void tearDown(void) { }

/*
 * CONVERT BCF:
 *  bcf_wahbm(in, out, bim, vid, tmp_dir, num_fields, num_records);
 *      convert_file_by_name_bcf_to_wahbm_bim(in,
 *                                            num_fields,
 *                                            num_records,
 *                                            wah_out,
 *                                            bim_out,
 *                                            vid_out,
 *                                            tmp_dir);
 *          init_bcf_file(in);
 *          push_bcf_gt_md(&q,
 *                         &bcf_f,
 *                         md_index,
 *                         num_inds,
 *                         num_vars,
 *                         gt_of_name,
 *                         md_of_name);
 *              bcf_get_genotypes(bcf_f->hdr,
 *                                bcf_f->line,
 *                                &gt_p,
 *                                &ntmp);
 *          sort_gt_md(&q,
 *                     md_index,
 *                     md_s_index,
 *                     num_inds,
 *                     num_vars,
 *                     gt_of_name,
 *                     gt_s_of_name,
 *                     md_of_name,
 *                     md_s_of_name,
 *                     vid_out);
 *          compress_md(&bcf_f,
 *                      md_of_name,
 *                      bim_out,
 *                      md_index,
 *                      num_vars);
 *          rotate_gt(num_inds,
 *                    num_vars,
 *                    gt_s_of_name,
 *                    gt_s_r_of_name);
 *          convert_file_by_name_ubin_to_wahbm(gt_s_r_of_name, wah_out);
 */

void test_init_bcf_file(void)
{
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);
    TEST_ASSERT_EQUAL(NUM_INDS, bcf_f.num_records);
    close_bcf_file(&bcf_f);
}

void test_push_bcf_gt_md(void)
{
    pri_queue q = priq_new(0);

    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);


    uint64_t *md_index = (uint64_t *) malloc(NUM_VARS * sizeof(uint64_t));

    char *gt_of_name = "tmp.gt_of_name";
    char *md_of_name = "tmp.md_of_name";

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   NUM_INDS,
                   NUM_VARS,
                   gt_of_name,
                   md_of_name);

    FILE *f = fopen(gt_of_name, "rb");

    


    //remove(gt_of_name);
    //remove(md_of_name);

    free(md_index);
}
