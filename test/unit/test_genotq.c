#include <stdio.h>
#include <stdlib.h>
#include "parse_q.h"
#include "bcf.h"
#include "off.h"
#include "wahbm.h"
#include "ubin.h"
#include "plt.h"
#include "genotq.h"
#include "unity.h"
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

char *BCF_FILE = "../data/10.1e4.var.bcf";
uint64_t BCF_OFFSETS[3] = {577,645,713};
uint64_t VCF_OFFSETS[3] = {670,737,804};

char *VCF_FILE = "../data/10.1e4.var.vcf.gz";
uint32_t NUM_INDS = 10;
uint32_t NUM_VARS = 43;

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
//{{{void test_init_bcf_file(void)
void test_init_bcf_file(void)
{
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);
    TEST_ASSERT_EQUAL(NUM_INDS, bcf_f.num_records);
    TEST_ASSERT_EQUAL(1, bcf_f.is_bcf);
    close_bcf_file(&bcf_f);

    struct bcf_file vcf_f = init_bcf_file(VCF_FILE);
    TEST_ASSERT_EQUAL(NUM_INDS, vcf_f.num_records);
    TEST_ASSERT_EQUAL(0, vcf_f.is_bcf);
    close_bcf_file(&vcf_f);
}
//}}}

//{{{void test_get_bcf_line(void)
void test_get_bcf_line(void)
{
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);

    TEST_ASSERT_EQUAL(BCF_OFFSETS[0], bcf_f.offset);

    TEST_ASSERT_EQUAL(0, get_bcf_line(&bcf_f));
    TEST_ASSERT_EQUAL(0, bcf_f.line->rid);
    TEST_ASSERT_EQUAL(0, bcf_f.line->pos);
    TEST_ASSERT_EQUAL(1, bcf_f.line->rlen);

    TEST_ASSERT_EQUAL(BCF_OFFSETS[1], bcf_f.offset);

    TEST_ASSERT_EQUAL(0, get_bcf_line(&bcf_f));
    TEST_ASSERT_EQUAL(0, bcf_f.line->rid);
    TEST_ASSERT_EQUAL(1, bcf_f.line->pos);
    TEST_ASSERT_EQUAL(1, bcf_f.line->rlen);

    TEST_ASSERT_EQUAL(BCF_OFFSETS[2], bcf_f.offset);

    close_bcf_file(&bcf_f);

    struct bcf_file vcf_f = init_bcf_file(VCF_FILE);

    //TEST_ASSERT_EQUAL(VCF_OFFSETS[0], vcf_f.offset);

    TEST_ASSERT_EQUAL(0, get_bcf_line(&vcf_f));
    TEST_ASSERT_EQUAL(0, vcf_f.line->rid);
    TEST_ASSERT_EQUAL(0, vcf_f.line->pos);
    TEST_ASSERT_EQUAL(1, vcf_f.line->rlen);

    //TEST_ASSERT_EQUAL(VCF_OFFSETS[1], vcf_f.offset);

    TEST_ASSERT_EQUAL(0, get_bcf_line(&vcf_f));
    TEST_ASSERT_EQUAL(0, vcf_f.line->rid);
    TEST_ASSERT_EQUAL(1, vcf_f.line->pos);
    TEST_ASSERT_EQUAL(1, vcf_f.line->rlen);

    //TEST_ASSERT_EQUAL(VCF_OFFSETS[2], vcf_f.offset);

    close_bcf_file(&vcf_f);
}
//}}}

//{{{void test_push_bcf_gt_offset(void)
void test_push_bcf_gt_offset(void)
{
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);
    pri_queue q_bcf = priq_new(0);

    push_bcf_gt_offset(&q_bcf,
                       &bcf_f,
                       NUM_INDS,
                       NUM_VARS,
                       ".tmp.gt_file",
                       ".tmp.offset_file",
                       "test_push_bcf_gt_offset");

    close_bcf_file(&bcf_f);

    struct off_file *of_bcf = open_off_file(".tmp.offset_file");

    TEST_ASSERT_EQUAL(BCF_OFFSETS[0], of_bcf->offsets[0]);
    TEST_ASSERT_EQUAL(BCF_OFFSETS[1], of_bcf->offsets[1]);
    TEST_ASSERT_EQUAL(BCF_OFFSETS[2], of_bcf->offsets[2]);

    destroy_off_file(of_bcf);

    struct bcf_file vcf_f = init_bcf_file(VCF_FILE);
    pri_queue q_vcf = priq_new(0);

    push_bcf_gt_offset(&q_vcf,
                       &vcf_f,
                       NUM_INDS,
                       NUM_VARS,
                       ".tmp.gt_file",
                       ".tmp.offset_file",
                       "test_push_bcf_gt_offset");

    close_bcf_file(&vcf_f);

    struct off_file *of_vcf = open_off_file(".tmp.offset_file");

    //TEST_ASSERT_EQUAL(VCF_OFFSETS[0], of_vcf->offsets[0]);
    //TEST_ASSERT_EQUAL(VCF_OFFSETS[1], of_vcf->offsets[1]);
    //TEST_ASSERT_EQUAL(VCF_OFFSETS[2], of_vcf->offsets[2]);

    destroy_off_file(of_vcf);
}
//}}}

//{{{void test_push_bcf_gt_md(void)
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
//}}}

//{{{void test_append_wahbm_to_wahbm_file(void)
void test_plt_line_to_packed_ints(void)
{

    uint32_t num_variants = 60;

    char *plt_1 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_2 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_3 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";

    char *plt_4 = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";


    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, num_variants, &ubin_1);
    TEST_ASSERT_EQUAL(1 + ((num_variants-1)/16), ubin_len_1);
    TEST_ASSERT_EQUAL( 0, ubin_1[0]);
    TEST_ASSERT_EQUAL( 0, ubin_1[1]);
    TEST_ASSERT_EQUAL( 0, ubin_1[2]);
    TEST_ASSERT_EQUAL( 0, ubin_1[3]);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, num_variants, &ubin_2);
    TEST_ASSERT_EQUAL(1 + ((num_variants-1)/16), ubin_len_2);
    // 01010101010101010101010101010101
    TEST_ASSERT_EQUAL( 1431655765, ubin_2[0]);
    // 01010101010101010101010101010101
    TEST_ASSERT_EQUAL( 1431655765, ubin_2[1]);
    // 01010101010101010101010101010101
    TEST_ASSERT_EQUAL( 1431655765, ubin_2[2]);
    // 01010101010101010101010100000000
    TEST_ASSERT_EQUAL( 1431655680, ubin_2[3]);

    uint32_t *ubin_3;
    uint32_t ubin_len_3 = plt_line_to_packed_ints(plt_3, num_variants, &ubin_3);
    TEST_ASSERT_EQUAL(1 + ((num_variants-1)/16), ubin_len_3);
    //01000100010001000100010001000100
    TEST_ASSERT_EQUAL( 1145324612, ubin_3[0]);
    //01000100010001000100010001000100
    TEST_ASSERT_EQUAL( 1145324612, ubin_3[1]);
    //01000100010001000100010001000100
    TEST_ASSERT_EQUAL( 1145324612, ubin_3[2]);
    //01000100010001000100010000000000
    TEST_ASSERT_EQUAL( 1145324544, ubin_3[3]);

    uint32_t *ubin_4;
    uint32_t ubin_len_4 = plt_line_to_packed_ints(plt_4, num_variants, &ubin_4);
    TEST_ASSERT_EQUAL(1 + ((num_variants-1)/16), ubin_len_4);
    //01000000000000000000000000000000
    TEST_ASSERT_EQUAL( 1073741824, ubin_4[0]);
    //00000000000000000000000000000100
    TEST_ASSERT_EQUAL( 4, ubin_4[1]);
    //00000000000000000000000000000000
    TEST_ASSERT_EQUAL( 0, ubin_4[2]);
    //00000000000000000000000000000000
    TEST_ASSERT_EQUAL( 0, ubin_4[3]);



    /*
    uint32_t *wah_1;
    uint32_t *wah_sizes_1;
    uint32_t wah_len_1 = ubin_to_bitmap_wah(ubin_1,
                                            ubin_len_1,
                                            60,
                                            &wah_1,
                                            &wah_sizes_1);
 
    */
}
//}}}

//{{{ void test_ubin_to_bitmap_wah(void)
void test_ubin_to_bitmap_wah(void)
{

    uint32_t num_variants = 60;

    char *plt_1 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_2 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_3 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, num_variants, &ubin_1);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, num_variants, &ubin_2);

    uint32_t *ubin_3;
    uint32_t ubin_len_3 = plt_line_to_packed_ints(plt_3, num_variants, &ubin_3);

    uint32_t *wah_1;
    uint32_t *wah_sizes_1;
    uint32_t wah_len_1 = ubin_to_bitmap_wah(ubin_1,
                                            ubin_len_1,
                                            60,
                                            &wah_1,
                                            &wah_sizes_1);
    TEST_ASSERT_EQUAL( 2, wah_sizes_1[0]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_1[1]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_1[2]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_1[3]);

    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111111"),
                       wah_1[0]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111100"),
                       wah_1[1]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[2]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[3]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[4]);

    /************************************************************************/
    uint32_t *wah_2;
    uint32_t *wah_sizes_2;
    uint32_t wah_len_2 = ubin_to_bitmap_wah(ubin_2,
                                            ubin_len_2,
                                            60,
                                            &wah_2,
                                            &wah_sizes_2);
    TEST_ASSERT_EQUAL( 1, wah_sizes_2[0]);
    TEST_ASSERT_EQUAL( 2, wah_sizes_2[1]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_2[2]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_2[3]);

    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_2[0]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111111"),
                       wah_2[1]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111100"),
                       wah_2[2]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_2[3]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_2[4]);

    /************************************************************************/

    uint32_t *wah_3;
    uint32_t *wah_sizes_3;
    uint32_t wah_len_3 = ubin_to_bitmap_wah(ubin_3,
                                            ubin_len_3,
                                            60,
                                            &wah_3,
                                            &wah_sizes_3);

    TEST_ASSERT_EQUAL( 2, wah_sizes_3[0]);
    TEST_ASSERT_EQUAL( 2, wah_sizes_3[1]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_3[2]);
    TEST_ASSERT_EQUAL( 1, wah_sizes_3[3]);

    TEST_ASSERT_EQUAL( bin_char_to_int("00101010101010101010101010101010"),
                       wah_3[0]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01010101010101010101010101010100"),
                       wah_3[1]);

    TEST_ASSERT_EQUAL( bin_char_to_int("01010101010101010101010101010101"),
                       wah_3[2]);
    TEST_ASSERT_EQUAL( bin_char_to_int("00101010101010101010101010101000"),
                       wah_3[3]);

    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_3[4]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_3[5]);
}
//}}}

//{{{ void test_convert_file_by_name_ubin_to_wahbm(void)
void test_convert_file_by_name_ubin_to_wahbm(void)
{
    char *ubin_file_name = "test_ubin";
    char *gqt_file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 60;
    uint32_t num_samples = 3;

    char *plt_0 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_1 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_2 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";

    uint32_t A[16] = {bin_char_to_int("01111111111111111111111111111111"),
                      bin_char_to_int("01111111111111111111111111111100"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("01111111111111111111111111111111"),
                      bin_char_to_int("01111111111111111111111111111100"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("00101010101010101010101010101010"),
                      bin_char_to_int("01010101010101010101010101010100"),
                      bin_char_to_int("01010101010101010101010101010101"),
                      bin_char_to_int("00101010101010101010101010101000"),
                      bin_char_to_int("10000000000000000000000000000010"),
                      bin_char_to_int("10000000000000000000000000000010")};

    uint32_t *ubin_0;
    uint32_t ubin_len_0 = plt_line_to_packed_ints(plt_0, 60, &ubin_0);

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, 60, &ubin_1);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, 60, &ubin_2);

    FILE *f = fopen(ubin_file_name, "wb");
    if (!f)
        err(EX_IOERR, "Error opening to \"%s\"", ubin_file_name); 

    if (fwrite(&num_variants, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", ubin_file_name); 

    if (fwrite(&num_samples, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", ubin_file_name); 

    if (fwrite(ubin_0, sizeof(uint32_t), ubin_len_0, f) != ubin_len_0)
        err(EX_IOERR, "Error writing to \"%s\"", ubin_file_name); 

    if (fwrite(ubin_1, sizeof(uint32_t), ubin_len_1, f) != ubin_len_1)
        err(EX_IOERR, "Error writing to \"%s\"", ubin_file_name); 

    if (fwrite(ubin_2, sizeof(uint32_t), ubin_len_2, f) != ubin_len_2)
        err(EX_IOERR, "Error writing to \"%s\"", ubin_file_name); 

    fclose(f);

    uint32_t r = convert_file_by_name_ubin_to_wahbm(ubin_file_name,
                                                    gqt_file_name,
                                                    full_cmd);

    struct wahbm_file *o = open_wahbm_file(gqt_file_name);

    uint32_t *wah_bitmap = NULL;;
    uint32_t wah_size = get_wahbm_bitmap(o, 0, 0, &wah_bitmap);

    TEST_ASSERT_EQUAL(2, wah_size);
    TEST_ASSERT_EQUAL(A[0], wah_bitmap[0]);
    TEST_ASSERT_EQUAL(A[1], wah_bitmap[1]);

    wah_size = get_wahbm_bitmap(o, 0, 1, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[2], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 0, 2, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[3], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 0, 3, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[4], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 1, 0, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[5], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 1, 1, &wah_bitmap);
    TEST_ASSERT_EQUAL(2, wah_size);
    TEST_ASSERT_EQUAL(A[6], wah_bitmap[0]);
    TEST_ASSERT_EQUAL(A[7], wah_bitmap[1]);

    wah_size = get_wahbm_bitmap(o, 1, 2, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[8], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 1, 3, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[9], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 2, 0, &wah_bitmap);
    TEST_ASSERT_EQUAL(2, wah_size);
    TEST_ASSERT_EQUAL(A[10], wah_bitmap[0]);
    TEST_ASSERT_EQUAL(A[11], wah_bitmap[1]);

    wah_size = get_wahbm_bitmap(o, 2, 1, &wah_bitmap);
    TEST_ASSERT_EQUAL(2, wah_size);
    TEST_ASSERT_EQUAL(A[12], wah_bitmap[0]);
    TEST_ASSERT_EQUAL(A[13], wah_bitmap[1]);

    wah_size = get_wahbm_bitmap(o, 2, 2, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[14], wah_bitmap[0]);

    wah_size = get_wahbm_bitmap(o, 2, 3, &wah_bitmap);
    TEST_ASSERT_EQUAL(1, wah_size);
    TEST_ASSERT_EQUAL(A[15], wah_bitmap[0]);

    free(wah_bitmap);
    destroy_wahbm_file(o);
}
//}}}
