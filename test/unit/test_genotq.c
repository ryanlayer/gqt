#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <sys/stat.h>
#include <sqlite3.h>

#include "unity.h"
#include "off.h"
#include "bcf.h"
#include "bim.h"
#include "bm.h"
#include "genotq.h"
#include "output_buffer.h"
#include "parse_q.h"
#include "ped.h"
#include "plt.h"
#include "pq.h"
#include "pthread_pool.h"
#include "quick_file.h"
#include "timer.h"
#include "ubin.h"
#include "variant_metadata.h"
#include "vid.h"
#include "wah.h"
#include "wahbm.h"
#include "wahbm_compressed_in_place.h"
#include "wahbm_in_place.h"

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

    


    remove(gt_of_name);
    remove(md_of_name);
    free(md_index);
}
//}}}

//{{{void test_uint32_t_to_bitmap(void)
void test_uint32_t_to_bitmap(void)
{
    uint32_t I[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10 };

    uint32_t bin_range_lo_I[3] = { 3, 6, 7 };
    uint32_t bin_range_hi_I[3] = { 5, 7, 9 };

    uint32_t *bm = NULL;

    uint32_t bit_array_lens = uint32_t_to_bitmap(I,
                                                 10,
                                                 bin_range_lo_I,
                                                 bin_range_hi_I,
                                                 3,
                                                 1,
                                                 1,
                                                 5,
                                                 &bm);
    TEST_ASSERT_EQUAL(1, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);

    // everything < 3: 0,1
    TEST_ASSERT_TRUE(bm[0] == (1 << (31 - 0)) + (1 << (31 - 1)));
    // everything [3,5): 2,3
    TEST_ASSERT_TRUE(bm[1] == (1 << (31 - 2)) + (1 << (31 - 3)));
    // everything [6,7): 5
    TEST_ASSERT_TRUE(bm[2] == (1 << (31 - 5)));
    // everything [7,9): 6,7
    TEST_ASSERT_TRUE(bm[3] == (1 << (31 - 6)) + (1 << (31 - 7)));
    // everything >=9: 8,9
    TEST_ASSERT_TRUE(bm[4] == (1 << (31 - 8)) + (1 << (31 - 9)));


    // test without upper and lower bins
    free(bm);
    bm = NULL;
    bit_array_lens = uint32_t_to_bitmap(I,
                                        10,
                                        bin_range_lo_I,
                                        bin_range_hi_I,
                                        3,
                                        0,
                                        0,
                                        5,
                                        &bm);
    TEST_ASSERT_EQUAL(1, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);

    // everything [3,5): 2,3
    TEST_ASSERT_TRUE(bm[0] == (1 << (31 - 2)) + (1 << (31 - 3)));
    // everything [6,7): 5
    TEST_ASSERT_TRUE(bm[1] == (1 << (31 - 5)));
    // everything [7,9): 6,7
    TEST_ASSERT_TRUE(bm[2] == (1 << (31 - 6)) + (1 << (31 - 7)));

    // test something larger
    free(bm);
    bm = NULL;
    uint32_t *J = (uint32_t *)malloc(100 * sizeof(uint32_t));

    uint32_t i,j;
    for (i = 0; i < 100; ++i) 
        J[i] = rand() % 100;

    uint32_t bin_range_lo_J[4] = {10, 30, 50, 80};
    uint32_t bin_range_hi_J[4] = {20, 40, 70, 90};

    bit_array_lens = uint32_t_to_bitmap(J,
                                        100,
                                        bin_range_lo_J,
                                        bin_range_hi_J,
                                        4,
                                        1,
                                        0,
                                        5,
                                        &bm);

    TEST_ASSERT_EQUAL(4, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);


    int bit_i = 0, int_i = 0, num_set = 0, bm_i = 0;
    uint32_t *target_bm;
    for (i = 0; i < 100; ++i) {

        num_set = 0;
        bm_i = 0;
        for (j = 0; j < 5; ++j) {
            target_bm = bm + j*bit_array_lens;
            if ( ((target_bm[int_i] >> (31 - bit_i)) & 1) == 1) {
                num_set += 1;
                bm_i = j;
            }
        }

        TEST_ASSERT_TRUE(num_set < 2);

        if ((num_set == 1) && (bm_i > 0)) {
            TEST_ASSERT_TRUE(J[i] >= bin_range_lo_J[bm_i-1]);
            TEST_ASSERT_TRUE(J[i] <  bin_range_hi_J[bm_i-1]);
        } else if ((num_set == 1) && (bm_i == 0)) {
            TEST_ASSERT_TRUE(J[i] < 10);
        }

        bit_i += 1;
        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    free(bm);
    free(J);
}
//}}}

//{{{void test_double_to_bitmap(void)
void test_double_to_bitmap(void)
{
    double I[10] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 , 1.0 };

    double bin_range_lo_I[3] = { 0.3, 0.6, 0.7 };
    double bin_range_hi_I[3] = { 0.5, 0.7, 0.9 };

    uint32_t *bm = NULL;

    uint32_t bit_array_lens = double_to_bitmap(I,
                                               10,
                                               bin_range_lo_I,
                                               bin_range_hi_I,
                                               3,
                                               1,
                                               1,
                                               5,
                                               &bm);
    TEST_ASSERT_EQUAL(1, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);

    // everything < 3: 0,1
    TEST_ASSERT_TRUE(bm[0] == (1 << (31 - 0)) + (1 << (31 - 1)));
    // everything [3,5): 2,3
    TEST_ASSERT_TRUE(bm[1] == (1 << (31 - 2)) + (1 << (31 - 3)));
    // everything [6,7): 5
    TEST_ASSERT_TRUE(bm[2] == (1 << (31 - 5)));
    // everything [7,9): 6,7
    TEST_ASSERT_TRUE(bm[3] == (1 << (31 - 6)) + (1 << (31 - 7)));
    // everything >=9: 8,9
    TEST_ASSERT_TRUE(bm[4] == (1 << (31 - 8)) + (1 << (31 - 9)));

    // test without upper and lower bins
    free(bm);
    bm = NULL;
    bit_array_lens = double_to_bitmap(I,
                                      10,
                                      bin_range_lo_I,
                                      bin_range_hi_I,
                                      3,
                                      0,
                                      0,
                                      5,
                                      &bm);
    TEST_ASSERT_EQUAL(1, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);

    // everything [3,5): 2,3
    TEST_ASSERT_TRUE(bm[0] == (1 << (31 - 2)) + (1 << (31 - 3)));
    // everything [6,7): 5
    TEST_ASSERT_TRUE(bm[1] == (1 << (31 - 5)));
    // everything [7,9): 6,7
    TEST_ASSERT_TRUE(bm[2] == (1 << (31 - 6)) + (1 << (31 - 7)));

    // test something larger
    free(bm);
    bm = NULL;
    double *J = (double *)malloc(100 * sizeof(double));

    uint32_t i,j;
    for (i = 0; i < 100; ++i) 
        J[i] = ((double)(rand() % 100))/(100.0);

    double bin_range_lo_J[4] = {0.10, 0.30, 0.50, 0.80};
    double bin_range_hi_J[4] = {0.20, 0.40, 0.70, 0.90};

    bit_array_lens = double_to_bitmap(J,
                                      100,
                                      bin_range_lo_J,
                                      bin_range_hi_J,
                                      4,
                                      1,
                                      0,
                                      5,
                                      &bm);

    TEST_ASSERT_EQUAL(4, bit_array_lens);
    TEST_ASSERT_TRUE(bm != NULL);


    int bit_i = 0, int_i = 0, num_set = 0, bm_i = 0;
    uint32_t *target_bm;
    for (i = 0; i < 100; ++i) {

        num_set = 0;
        bm_i = 0;
        for (j = 0; j < 5; ++j) {
            target_bm = bm + j*bit_array_lens;
            if ( ((target_bm[int_i] >> (31 - bit_i)) & 1) == 1) {
                num_set += 1;
                bm_i = j;
            }
        }

        TEST_ASSERT_TRUE(num_set < 2);

        if ((num_set == 1) && (bm_i > 0)) {
            TEST_ASSERT_TRUE(J[i] >= bin_range_lo_J[bm_i-1]);
            TEST_ASSERT_TRUE(J[i] <  bin_range_hi_J[bm_i-1]);
        } else if ((num_set == 1) && (bm_i == 0)) {
            TEST_ASSERT_TRUE(J[i] < 10);
        }

        bit_i += 1;
        if (bit_i == 32) {
            bit_i = 0;
            int_i += 1;
        }
    }

    free(bm);
    free(J);
}
//}}}

//{{{ void test_ints_to_wah(void)
void test_ints_to_wah(void)
{
    /*
     * |-32---------------------------|
     * 10000000000000000000000000000000 -> 2147483648
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000011 -> 3
     * 00000000000000000000000000000001 -> 1
     *
     * Regroup into 31-bit groups
     *
     * |-31--------------------------|
     * 1000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0011000000000000000000000000000
     * 00001
     *
     * Pad them out
     *
     * |-31--------------------------|
     * 1000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0011000000000000000000000000000
     * 0000100000000000000000000000000
     *
     * First one will be represented as a litteral
     * 01000000000000000000000000000000
     * The next three will compress to a fill of 93 zeros (3 words)
     * 10000000000000000000000000000011
     * Next two will also be litteral
     * 00011000000000000000000000000000
     * 00000100000000000000000000000000
     *
     * 01000000000000000000000000000000 -> 0x40000000
     * 10000000000000000000000000000011 -> 0x80000003
     * 00011000000000000000000000000000 -> 0x18000000
     * 00000100000000000000000000000000 -> 0x04000000
     *
     */
    unsigned int I[5] = {2147483648,0,0,3,1};
    //unsigned int A[4] = {1073741824, 2147483741, 402653184, 67108864};
    unsigned int A[4] = {0x40000000, 0x80000003, 0x18000000, 0x04000000};

    unsigned int *O;
    unsigned int wah_size = ints_to_wah(I,5,160,&O);

    TEST_ASSERT_EQUAL(4, wah_size);

    int i;
    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A[i]);

    /*
     * Move to 32 bits
     *
     * |-32---------------------------|
     * 01101110110001101110111010101010 -> 1858530986
     * 10101111000101011010101010101001 -> 2937432745
     * 01011110001010011011011001010101 -> 1579791957
     * 11010101110101010010111011010101 -> 3587518165
     * 01011010100000001010100001010101 -> 1518381141
     * 01010101010000000000000000000000 -> 1430257664
     *
     * |-31--------------------------|
     * 0110111011000110111011101010101
     * 0101011110001010110101010101010
     * 0101011110001010011011011001010
     * 1011101010111010101001011101101
     * 0101010110101000000010101000010
     * 1010101010101010000000000000000
     * 000000
 
     * 31 with padding
     * |-31--------------------------|
     * 00110111011000110111011101010101 -> 929265493
     * 00101011110001010110101010101010 -> 734358186
     * 00101011110001010011011011001010 -> 734344906
     * 01011101010111010101001011101101 -> 1566397165
     * 00101010110101000000010101000010 -> 718538050
     * 01010101010101010000000000000000 -> 1431633920
     * 00000000000000000000000000000000 -> 0
     */

    unsigned int I2[6] = {1858530986,
                          2937432745,
                          1579791957,
                          3587518165,
                          1518381141,
                          1430257664};
    unsigned int A2[7] = {929265493,
                          734358186,
                          734344906,
                          1566397165,
                          718538050,
                          1431633920,
                          0};
    
    free(O);

    wah_size = ints_to_wah(I2,6,192,&O);

    TEST_ASSERT_EQUAL(7, wah_size);

    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A2[i]);

    /*
     * |-32---------------------------|
     * 10000000000000000000000000000000 -> 
     * 11111111111111111111111111111111 -> 
     * 11111111111111111111111111111111 -> 
     * 11111111111111111111111111111111 -> 
     * 00000000000000000000000000000011 -> 
     * 00000000000000000000000000000001 ->
     *
     * Regroup into 31-bit groups
     * |-31--------------------------|
     * 1000000000000000000000000000000
     * 0111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1111000000000000000000000000000
     * 0001100000000000000000000000000
     * 000001
     *
     *
     * 01000000000000000000000000000000
     * 00111111111111111111111111111111
     * 11000000000000000000000000000010
     * 01111000000000000000000000000000
     * 00001100000000000000000000000000
     * 00000010000000000000000000000000
     */
    unsigned int I3[6] = {
            bin_char_to_int("10000000000000000000000000000000"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("00000000000000000000000000000011"),
            bin_char_to_int("00000000000000000000000000000001")
        };
 
    unsigned int A3[6] = {
            bin_char_to_int("01000000000000000000000000000000"),
            bin_char_to_int("00111111111111111111111111111111"),
            bin_char_to_int("11000000000000000000000000000010"),
            bin_char_to_int("01111000000000000000000000000000"),
            bin_char_to_int("00001100000000000000000000000000"),
            bin_char_to_int("00000010000000000000000000000000")
        };

    free(O);
    wah_size = ints_to_wah(I3,6,192,&O);

    TEST_ASSERT_EQUAL(6, wah_size);

    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A3[i]);
}
//}}}

//{{{void test_uint32_to_wah_bitmap(void)
void test_uint32_to_wah_bitmap(void)
{
    uint32_t *J = (uint32_t *)malloc(100 * sizeof(uint32_t));

    uint32_t i,j;
    for (i = 0; i < 100; ++i) 
        J[i] = 20 + (rand() % 50);

    uint32_t bin_range_lo_J[4] = {10, 30, 50, 80};
    uint32_t bin_range_hi_J[4] = {20, 40, 70, 90};

    uint32_t **ws = NULL;
    uint32_t *w_lens = NULL;

    uint32_t bit_array_lens = uint32_t_to_wah_bitmap(J,
                                                     100,
                                                     bin_range_lo_J,
                                                     bin_range_hi_J,
                                                     4,
                                                     1,
                                                     0,
                                                     5,
                                                     &ws,
                                                     &w_lens);

    uint32_t lens[5] = {0,0,0,0,0};
    for (i = 0; i < 5; ++i) {
        //printf("%u\n", w_lens[i]);
        for (j = 0; j < w_lens[i]; ++j) {
            if ( (ws[i][j] >> 31) == 1 )
                lens[i] += ((ws[i][j] << 1) >> 1);
            else 
                lens[i] += 1;
        }
    }

    for (i = 0; i < 5; ++i) 
        TEST_ASSERT_EQUAL(4, lens[i]);
}
//}}}

//{{{void test_int_to_wah_bitmap(void)
void test_int_to_wah_bitmap(void)
{
    int *J = (int *)malloc(100 * sizeof(int));

    uint32_t i,j;
    for (i = 0; i < 100; ++i) 
        J[i] = (rand() % 2);

    int bin_range_lo_J[4] = {0, 1};
    int bin_range_hi_J[4] = {1, 2};

    uint32_t **ws = NULL;
    uint32_t *w_lens = NULL;

    uint32_t bit_array_lens = int_to_wah_bitmap(J,
                                                100,
                                                bin_range_lo_J,
                                                bin_range_hi_J,
                                                2,
                                                0,
                                                0,
                                                2,
                                                &ws,
                                                &w_lens);

    uint32_t lens[2] = {0,0};
    for (i = 0; i < 2; ++i) {
        //printf("%u\n", w_lens[i]);
        for (j = 0; j < w_lens[i]; ++j) {
            if ( (ws[i][j] >> 31) == 1 )
                lens[i] += ((ws[i][j] << 1) >> 1);
            else 
                lens[i] += 1;
        }
    }

    for (i = 0; i < 2; ++i) 
        TEST_ASSERT_EQUAL(4, lens[i]);
}
//}}}

//{{{void test_float_to_wah_bitmap(void)
void test_float_to_wah_bitmap(void)
{
    float *J = (float *)malloc(100 * sizeof(float));

    uint32_t i,j;
    for (i = 0; i < 100; ++i) 
        J[i] = ((float)(rand()%100))/100.0;

    float bin_range_lo_J[4] = {0.1, 0.3, 0.5, 0.8};
    float bin_range_hi_J[4] = {0.3, 0.5, 0.8, 0.9};

    uint32_t **ws = NULL;
    uint32_t *w_lens = NULL;

    uint32_t bit_array_lens = float_to_wah_bitmap(J,
                                                  100,
                                                  bin_range_lo_J,
                                                  bin_range_hi_J,
                                                  4,
                                                  1,
                                                  1,
                                                  6,
                                                  &ws,
                                                  &w_lens);

    uint32_t bit_offset = 0;
    for (i = 0; i < w_lens[0]; ++i) {
        if (ws[0][i]>> 31 == 1) {
            bit_offset += 31*((ws[0][i]<<1)>>1);
        } else {
            uint32_t j;
            for (j = 0; j < 31; j++) {
                if ( ((ws[0][i]>>(30-j))&1) == 1) {
                    TEST_ASSERT_TRUE(J[j+bit_offset] < bin_range_lo_J[0]);
                }
            }
            bit_offset += 31;
        }
    }
}
//}}}

//{{{ void test_get_variant_metadata(void)
void test_get_variant_metadata(void)
{
    uint32_t r;
    
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);

    int missing = -1;
    int type;
    r = get_variant_metadata(&bcf_f,
                             NUM_VARS,
                             "SHOULD_FAIL",
                             NULL,
                             &missing,
                             &type);

    TEST_ASSERT_EQUAL(1, r);

    close_bcf_file(&bcf_f);

    bcf_f = init_bcf_file(BCF_FILE);

    float A_LOF[43] = {1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,
                     0,0,0,0,0,0,0,0,1};
    float *LOF;
    missing = 0;
    r = get_variant_metadata(&bcf_f,
                             NUM_VARS,
                             "LOF",
                             (void **)&LOF,
                             &missing,
                             &type);

    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(0, type);

    uint32_t i;
    for (i = 0; i < NUM_VARS; ++i) 
        TEST_ASSERT_EQUAL(A_LOF[i], LOF[i]);

    close_bcf_file(&bcf_f);



    bcf_f = init_bcf_file(BCF_FILE);

    float A_AC[43] = {6,1,3,2,2,3,4,13,13,3,12,12,12,2,12,3,8,
                    2,1,12,2,12,2,10,3,3,3,1,3,17,3,2,17,13,
                    4,3,4,3,3,2,3,4,2};
    float *AC;
    r = get_variant_metadata(&bcf_f,
                             NUM_VARS,
                             "AC",
                             (void **)&AC,
                             &missing,
                             &type);

    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(1, type);

    for (i = 0; i < NUM_VARS; ++i) 
        TEST_ASSERT_EQUAL(A_AC[i], AC[i]);

    close_bcf_file(&bcf_f);

    bcf_f = init_bcf_file(BCF_FILE);

    float A_AF[43] = {0.300000,0.050000,0.150000,0.100000,0.100000,
                      0.150000,0.200000,0.650000,0.650000,0.150000,
                      0.600000,0.600000,0.600000,0.100000,0.600000,
                      0.150000,0.400000,0.100000,0.050000,0.600000,
                      0.100000,0.600000,0.100000,0.500000,0.150000,
                      0.150000,0.150000,0.050000,0.150000,0.850000,
                      0.150000,0.100000,0.850000,0.650000,0.200000,
                      0.150000,0.200000,0.150000,0.150000,0.100000,
                      0.150000,0.200000,0.100000};
    float *AF;
    r = get_variant_metadata(&bcf_f,
                             NUM_VARS,
                             "AF",
                             (void **)&AF,
                             &missing,
                             &type);

    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(2, type);

    for (i = 0; i < NUM_VARS; ++i)
        TEST_ASSERT_EQUAL(A_AF[i], AF[i]);

    close_bcf_file(&bcf_f);
}
//}}}

//{{{ void test_int_equal_freq_binning(void)
void test_int_equal_freq_binning(void)
{
    uint32_t i;

    int A[43] = {6,1,3,2,2,3,4,13,13,3,12,12,12,2,12,3,8,
                 2,1,12,2,12,2,10,3,3,3,1,3,17,3,2,17,13,
                 4,3,4,3,3,2,3,4,2};

    int A_lo[5] = {1,2,4,6,12};
    int A_hi[5] = {2,4,6,12,14};

    int to_test = 20, num_bins = 5;

    int *bin_range_lo = NULL, *bin_range_hi = NULL;

    uint32_t actual_num_bins = int_equal_freq_binning(A,
                                                      to_test,
                                                      num_bins,
                                                      &bin_range_lo,
                                                      &bin_range_hi);

    TEST_ASSERT_EQUAL(num_bins, actual_num_bins);
    for (i = 0; i < actual_num_bins; ++i) {
        TEST_ASSERT_EQUAL(A_lo[i], bin_range_lo[i]);
        TEST_ASSERT_EQUAL(A_hi[i], bin_range_hi[i]);
    }


    int B[24] = {1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0};
    int B_lo[2] = {0,1};
    int B_hi[2] = {1,2};

    to_test = 20, num_bins = 5;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL, bin_range_hi = NULL;

    actual_num_bins = int_equal_freq_binning(B,
                                             to_test,
                                             num_bins,
                                             &bin_range_lo,
                                             &bin_range_hi);

    TEST_ASSERT_EQUAL(2, actual_num_bins);
    for (i = 0; i < actual_num_bins; ++i) {
        TEST_ASSERT_EQUAL(B_lo[i], bin_range_lo[i]);
        TEST_ASSERT_EQUAL(B_hi[i], bin_range_hi[i]);
    }
}
//}}}

//{{{ void test_float_equal_freq_binning(void)
void test_floatint_equal_freq_binning(void)
{
    uint32_t i;

    float A[43] = {6,1,3,2,2,3,4,13,13,3,12,12,12,2,12,3,8,
                 2,1,12,2,12,2,10,3,3,3,1,3,17,3,2,17,13,
                 4,3,4,3,3,2,3,4,2};

    float A_lo[5] = {1,2,4,6,12};
    float A_hi[5] = {2,4,6,12,14};

    int to_test = 20, num_bins = 5;

    float *bin_range_lo = NULL, *bin_range_hi = NULL;

    uint32_t actual_num_bins = float_equal_freq_binning(A,
                                                        to_test,
                                                        num_bins,
                                                        &bin_range_lo,
                                                        &bin_range_hi);

    TEST_ASSERT_EQUAL(num_bins, actual_num_bins);
    for (i = 0; i < actual_num_bins; ++i) {
        TEST_ASSERT_EQUAL(A_lo[i], bin_range_lo[i]);
        TEST_ASSERT_EQUAL(A_hi[i], bin_range_hi[i]);
    }


    float B[24] = {1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0};
    float B_lo[2] = {0,1};
    float B_hi[2] = {1,2};

    to_test = 20, num_bins = 5;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL, bin_range_hi = NULL;

    actual_num_bins = float_equal_freq_binning(B,
                                               to_test,
                                              num_bins,
                                              &bin_range_lo,
                                              &bin_range_hi);

    TEST_ASSERT_EQUAL(2, actual_num_bins);
    for (i = 0; i < actual_num_bins; ++i) {
        TEST_ASSERT_EQUAL(B_lo[i], bin_range_lo[i]);
        TEST_ASSERT_EQUAL(B_hi[i], bin_range_hi[i]);
    }
}
//}}}

//{{{ void test_register_variant_metadata_index(void)
void test_register_variant_metadata_index(void)
{
    char *test_db_name = ".tmp.test.db";
    remove(test_db_name);

    int row_id = register_variant_metadata_index("test_bcf",
                                                 test_db_name,
                                                 "FIRST",
                                                 10,
                                                 1,
                                                 5);

    TEST_ASSERT_EQUAL(1, row_id);

    // make sure the file was craeted
    struct stat buffer;
    int r = stat(test_db_name, &buffer);
    TEST_ASSERT_EQUAL(0, r);

    // should not insert twice
    row_id = register_variant_metadata_index("test_bcf",
                                             test_db_name,
                                             "FIRST",
                                             10,
                                             1,
                                             5);
    TEST_ASSERT_EQUAL(-1, row_id);

    // a second one should have the next id
    row_id = register_variant_metadata_index("test_bcf",
                                             test_db_name,
                                             "SECOND",
                                             10,
                                             1,
                                             5);

    TEST_ASSERT_EQUAL(2, row_id);

    remove(test_db_name);
}
//}}}

//{{{void test_add_variant_metadata_float_bins(void)
void test_add_variant_metadata_float_bins(void)
{
    /*
    uint32_t add_variant_metadata_float_bins(char *db_file_name,
                                             int rowid,
                                             float *float_bin_range_lo,
                                             float *float_bin_range_hi,
                                             int actual_num_bins,
                                             int less_than_bin,
                                             int greater_than_bin)
    */

    char *test_db_name = ".tmp.test.db";
    remove(test_db_name);

    int row_id = register_variant_metadata_index("test_bcf",
                                                 test_db_name,
                                                 "FIRST",
                                                 10,
                                                 1,
                                                 3);

    float float_bin_range_lo[3] = { 0.0, 1.0, 2.0};
    float float_bin_range_hi[3] = { 1.0, 2.0, 3.0};

    int rows_added = add_variant_metadata_float_bins(test_db_name,
                                                     row_id,
                                                     float_bin_range_lo,
                                                     float_bin_range_hi,
                                                     3,
                                                     0,
                                                     0);

    // since there is no greater than or less than bins, only 3 rows should
    // have been added
    TEST_ASSERT_EQUAL(3, rows_added);

    // check to see if 3 rows where inserted
    sqlite3 *db;
    char *err_msg;
    int rc = sqlite3_open(test_db_name, &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        exit(1);
    }

    char *q_check;
    int r = asprintf(&q_check,
                     "SELECT * FROM variant_md_bins WHERE md_rowid = %d",
                     row_id);

    uint32_t num_rows = 0;
    rc = sqlite3_exec(db, q_check, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(3, num_rows);

    // check to see if no rows with hi or lo bins set
    r = asprintf(&q_check,
                 "SELECT * FROM variant_md_bins WHERE "
                 "(md_rowid == %d) AND "
                 "(lo_bin == 1 OR hi_bin == 1)",
                 row_id);

    num_rows = 0;
    rc = sqlite3_exec(db, q_check, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(0, num_rows);

    // second test has hi and lo bins
    row_id = register_variant_metadata_index("test_bcf",
                                             test_db_name,
                                             "SECOND",
                                             10,
                                             1,
                                             5);

    rows_added = add_variant_metadata_float_bins(test_db_name,
                                                 row_id,
                                                 float_bin_range_lo,
                                                 float_bin_range_hi,
                                                 3,
                                                 1,
                                                 1);

    // since greater than or less than bins are set, 5 rows should have been
    // added
    TEST_ASSERT_EQUAL(5, rows_added);

    // check to see if 5 rows where inserted
    r = asprintf(&q_check,
                 "SELECT * FROM variant_md_bins WHERE md_rowid = %d",
                 row_id);

    num_rows = 0;
    rc = sqlite3_exec(db, q_check, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(5, num_rows);

    // check to that 2 rows are found
    r = asprintf(&q_check,
                 "SELECT * FROM variant_md_bins WHERE "
                 "(md_rowid == %d) AND "
                 "(lo_bin == 1 OR hi_bin == 1)",
                 row_id);

    num_rows = 0;
    rc = sqlite3_exec(db, q_check, count_callback, &num_rows, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(2, num_rows);

    // make sure lo bin id is 0 
    r = asprintf(&q_check,
                 "SELECT bin FROM variant_md_bins WHERE lo_bin == 1;");

    uint32_t val = 0;
    rc = sqlite3_exec(db, q_check, get_single_int_callback, &val, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(0, val);

   
    // make sure hi bin id is row_added-1
    r = asprintf(&q_check,
                 "SELECT bin FROM variant_md_bins WHERE hi_bin==1;");

    val = 0;
    rc = sqlite3_exec(db, q_check, get_single_int_callback, &val, &err_msg);
    if( rc != SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        exit(1);
    }
    TEST_ASSERT_EQUAL(rows_added - 1, val);


    sqlite3_close(db);
    remove(test_db_name);
}
//}}}

//{{{void test_index_variant_metadata(void)
void test_index_variant_metadata(void)
{
    char *gt_file_name = ".tmp.gqt",
         *bim_file_name = ".tmp.bim",
         *vid_file_name = ".tmp.vid",
         *tmp_dir = ".",
         *variant_index_file_name = ".tmp.gqtv",
         *variant_db_name = ".tmp.gqtv.db";

    char *full_cmd = "test_index_variant_metadata";

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);

    uint32_t r = convert_file_by_name_bcf_to_wahbm_bim(BCF_FILE,
                                                       NUM_INDS,
                                                       NUM_VARS,
                                                       gt_file_name,
                                                       bim_file_name,
                                                       vid_file_name,
                                                       tmp_dir,
                                                       full_cmd);


    char *LOF_field = "LOF";
    uint32_t num_to_test = 20, num_bins = 50;
    void *bin_range_lo = NULL, *bin_range_hi = NULL;
    int less_than_bin = -1, greater_than_bin = -1;

    int type = index_variant_metadata(BCF_FILE,
                                      vid_file_name,
                                      variant_db_name,
                                      variant_index_file_name,
                                      NUM_VARS,
                                      LOF_field,
                                      num_to_test,
                                      &num_bins,
                                      &bin_range_lo,
                                      &bin_range_hi,
                                      &less_than_bin,
                                      &greater_than_bin);

    TEST_ASSERT_EQUAL(0, type);
    TEST_ASSERT_EQUAL(2, num_bins);
    TEST_ASSERT_EQUAL(0, less_than_bin);
    TEST_ASSERT_EQUAL(0, greater_than_bin);
    TEST_ASSERT_EQUAL(0, ((float *)bin_range_lo)[0]);
    TEST_ASSERT_EQUAL(1, ((float *)bin_range_hi)[0]);
    TEST_ASSERT_EQUAL(1, ((float *)bin_range_lo)[1]);
    TEST_ASSERT_EQUAL(2, ((float *)bin_range_hi)[1]);

    char *AF_field = "AF";
    num_to_test = 20;
    num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AF_field,
                                  num_to_test,
                                  &num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    TEST_ASSERT_EQUAL(2, type);
    TEST_ASSERT_EQUAL(10, num_bins);
    TEST_ASSERT_EQUAL(1, less_than_bin);
    TEST_ASSERT_EQUAL(1, greater_than_bin);

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);


}
//}}}

//{{{ void test_get_variant_metadata_bin_info(void)
void test_get_variant_metadata_bin_info(void)
{
    char *gt_file_name = ".test_get_variant_metadata_bin_info.gqt",
         *bim_file_name = ".test_get_variant_metadata_bin_info.bim",
         *vid_file_name = ".test_get_variant_metadata_bin_info.vid",
         *tmp_dir = ".",
         *variant_index_file_name = ".test_get_variant_metadata_bin_info.gqtv",
         *variant_db_name = ".test_get_variant_metadata_bin_info.gqtv.db";

    char *full_cmd = "test_get_variant_metadata_bin_info";

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);

    uint32_t r = convert_file_by_name_bcf_to_wahbm_bim(BCF_FILE,
                                                       NUM_INDS,
                                                       NUM_VARS,
                                                       gt_file_name,
                                                       bim_file_name,
                                                       vid_file_name,
                                                       tmp_dir,
                                                       full_cmd);


    char *LOF_field = "LOF";
    uint32_t num_to_test = 20, lof_num_bins = 50;
    void *bin_range_lo = NULL, *bin_range_hi = NULL;
    int less_than_bin = -1, greater_than_bin = -1;

    int type = index_variant_metadata(BCF_FILE,
                                      vid_file_name,
                                      variant_db_name,
                                      variant_index_file_name,
                                      NUM_VARS,
                                      LOF_field,
                                      num_to_test,
                                      &lof_num_bins,
                                      &bin_range_lo,
                                      &bin_range_hi,
                                      &less_than_bin,
                                      &greater_than_bin);

    char *AF_field = "AF";
    num_to_test = 20;
    uint32_t af_num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AF_field,
                                  num_to_test,
                                  &af_num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    char *source_file;
    uint64_t offset;
    uint32_t num_bins;
    int md_type;
    int rowid = get_variant_metadata_bin_info(variant_db_name,
                                              "LOF",
                                              &offset,
                                              &num_bins,
                                              &md_type,
                                              &source_file);

    TEST_ASSERT_EQUAL(1,rowid);
    TEST_ASSERT_EQUAL(lof_num_bins,num_bins);

    rowid = get_variant_metadata_bin_info(variant_db_name,
                                          "AF",
                                          &offset,
                                          &num_bins,
                                          &md_type,
                                          &source_file);

    TEST_ASSERT_EQUAL(2,rowid);
    TEST_ASSERT_EQUAL(af_num_bins, num_bins);

    rowid = get_variant_metadata_bin_info(variant_db_name,
                                          "SHOULD_FAIL",
                                          &offset,
                                          &num_bins,
                                          &md_type,
                                          &source_file);

    TEST_ASSERT_EQUAL(-1,rowid);

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);
}
//}}}

//{{{ void test_get_query_bins(void)
void test_get_query_bins(void)
{
    char *gt_file_name = ".test_get_query_bins.gqt",
         *bim_file_name = ".test_get_query_bins.bim",
         *vid_file_name = ".test_get_query_bins.vid",
         *tmp_dir = ".",
         *variant_index_file_name = ".test_get_query_bins.gqtv",
         *variant_db_name = ".test_get_query_bins.gqtv.db";

    char *full_cmd = "test_get_query_bins";

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);

    //{{{ Index LOF, AF, and AC
    uint32_t r = convert_file_by_name_bcf_to_wahbm_bim(BCF_FILE,
                                                       NUM_INDS,
                                                       NUM_VARS,
                                                       gt_file_name,
                                                       bim_file_name,
                                                       vid_file_name,
                                                       tmp_dir,
                                                       full_cmd);


    char *LOF_field = "LOF";
    uint32_t num_to_test = 20, lof_num_bins = 50;
    void *bin_range_lo = NULL, *bin_range_hi = NULL;
    int less_than_bin = -1, greater_than_bin = -1;

    int type = index_variant_metadata(BCF_FILE,
                                      vid_file_name,
                                      variant_db_name,
                                      variant_index_file_name,
                                      NUM_VARS,
                                      LOF_field,
                                      num_to_test,
                                      &lof_num_bins,
                                      &bin_range_lo,
                                      &bin_range_hi,
                                      &less_than_bin,
                                      &greater_than_bin);

    char *AF_field = "AF";
    num_to_test = 20;
    uint32_t af_num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AF_field,
                                  num_to_test,
                                  &af_num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    char *AC_field = "AC";
    num_to_test = 20;
    uint32_t ac_num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AC_field,
                                  num_to_test,
                                  &ac_num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    char *source_file;
    uint64_t offset;
    uint32_t num_bins;
    int af_md_type;
    int af_rowid = get_variant_metadata_bin_info(variant_db_name,
                                                 "AF",
                                                 &offset,
                                                 &num_bins,
                                                 &af_md_type,
                                                 &source_file);
    //}}}

    //{{{ TEST AF (float)
    float lo = 0.1, hi = 0.5;
    uint32_t *bins = NULL;
    uint32_t i;
    // test ==
    int num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                         "AF",
                                                         af_rowid,
                                                         num_bins,
                                                         af_md_type,
                                                         1,
                                                         0.1,
                                                         0.5,
                                                         &bins);
                   
    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(2,bins[0]);

    free(bins);
    bins = NULL;

    //{{{ < and <= test
    // test <
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     2,
                                                     0.1,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(2,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);


    free(bins);
    bins = NULL;

    // test <
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     2,
                                                     0.11,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(3,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);


    free(bins);
    bins = NULL;
    // test < when hitting lo bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     2,
                                                     0.001,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0,bins[0]);

    free(bins);
    bins = NULL;

    // test < when hitting hi bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     2,
                                                     100.001,
                                                     0.5,
                                                     &bins);
    TEST_ASSERT_EQUAL(num_bins,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test <=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     4,
                                                     0.1,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(3,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test <=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     4,
                                                     0.11,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(3,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test <= when hiting lo bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     4,
                                                     0.001,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0,bins[0]);

    free(bins);
    bins = NULL;

    // test <= when hiting hi bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     4,
                                                     100.001,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(num_bins,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;
    //}}}

    //{{{ > and >= test
    // test >
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     8,
                                                     0.1,
                                                     0.1,
                                                     &bins);

    TEST_ASSERT_EQUAL(8,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(2+i, bins[i]);


    free(bins);
    bins = NULL;

    // test >
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     8,
                                                     0.1,
                                                     0.11,
                                                     &bins);

    TEST_ASSERT_EQUAL(8,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(2+i, bins[i]);


    free(bins);
    bins = NULL;


    // test > when hitting lo bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     8,
                                                     0.001,
                                                     0.001,
                                                     &bins);

    TEST_ASSERT_EQUAL(num_bins,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test > when hitting hi bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     8,
                                                     100.001,
                                                     100.5,
                                                     &bins);
    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(9, bins[0]);

    free(bins);
    bins = NULL;

    // test >=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     16,
                                                     0.1,
                                                     0.1,
                                                     &bins);

    TEST_ASSERT_EQUAL(8,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(2+i, bins[i]);


    free(bins);
    bins = NULL;

    // test >=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     16,
                                                     0.1,
                                                     0.11,
                                                     &bins);

    TEST_ASSERT_EQUAL(8,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(2+i, bins[i]);


    free(bins);
    bins = NULL;


    // test >= when hitting lo bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     16,
                                                     0.001,
                                                     0.001,
                                                     &bins);

    TEST_ASSERT_EQUAL(num_bins,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test >= when hitting hi bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     16,
                                                     100.001,
                                                     100.5,
                                                     &bins);
    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(9, bins[0]);

    free(bins);
    bins = NULL;



#if 0
    // test >=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     rowid,
                                                     num_bins,
                                                     16,
                                                     0.1,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(3,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;

    // test <= when hiting lo bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     rowid,
                                                     num_bins,
                                                     4,
                                                     0.001,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0,bins[0]);

    free(bins);
    bins = NULL;

    // test <= when hiting hi bin
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     rowid,
                                                     num_bins,
                                                     4,
                                                     100.001,
                                                     0.5,
                                                     &bins);

    TEST_ASSERT_EQUAL(num_bins,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_EQUAL(i, bins[i]);

    free(bins);
    bins = NULL;
#endif
    //}}}

    // test !=
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AF",
                                                     af_rowid,
                                                     num_bins,
                                                     af_md_type,
                                                     32,
                                                     lo,
                                                     hi,
                                                     &bins);

    TEST_ASSERT_EQUAL(num_bins - 1,num_query_bins);
    for (i = 0; i < num_query_bins; ++i)
        TEST_ASSERT_FALSE(bins[i] == 2);

    free(bins);
    bins = NULL;

    //}}}
   
    int ac_md_type;
    int ac_rowid = get_variant_metadata_bin_info(variant_db_name,
                                                 "AC",
                                                 &offset,
                                                 &num_bins,
                                                 &ac_md_type,
                                                 &source_file);

    //{{{ TEST AC
    /*
     * select * from variant_md_bins where md_rowid==3 order by bin;
     *  3	0	0.0	1.0	1	0
     *  3	1	1.0	2.0	0	0
     *  3	2	2.0	3.0	0	0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     *  3	7	12.0	13.0	0	0
     *  3	8	13.0	14.0	0	0
     *  3	9	14.0	0.0	0	1
     */

    /*
     * == 2
     * 3       2       2.0     3.0     0       0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     1,
                                                     2.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(2,bins[0]);

    free(bins);
    bins = NULL;

    /* < 2
     *  3	0	0.0	1.0	1	0
     *  3	1	1.0	2.0	0	0
     */ 
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     2,
                                                     2.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(2,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);
    TEST_ASSERT_EQUAL(1, bins[1]);

    free(bins);
    bins = NULL;

    /*
     * < when hitting lo bin
     *  3	0	0.0	1.0	1	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     2,
                                                     -2.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);

    free(bins);
    bins = NULL;

    /*
     * < when hitting hi bin
     *  3	0	0.0	1.0	1	0
     *  3	1	1.0	2.0	0	0
     *  3	2	2.0	3.0	0	0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     *  3	7	12.0	13.0	0	0
     *  3	8	13.0	14.0	0	0
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     2,
                                                     20.0,
                                                     0.0,
                                                     &bins);
    TEST_ASSERT_EQUAL(10,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);
    TEST_ASSERT_EQUAL(1, bins[1]);
    TEST_ASSERT_EQUAL(2, bins[2]);
    TEST_ASSERT_EQUAL(3, bins[3]);
    TEST_ASSERT_EQUAL(4, bins[4]);
    TEST_ASSERT_EQUAL(5, bins[5]);
    TEST_ASSERT_EQUAL(6, bins[6]);
    TEST_ASSERT_EQUAL(7, bins[7]);
    TEST_ASSERT_EQUAL(8, bins[8]);
    TEST_ASSERT_EQUAL(9, bins[9]);

    free(bins);
    bins = NULL;

    /*
     * <= 2.0
     *  3	0	0.0	1.0	1	0
     *  3	1	1.0	2.0	0	0
     *  3	2	2.0	3.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     4,
                                                     2.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(3,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);
    TEST_ASSERT_EQUAL(1, bins[1]);
    TEST_ASSERT_EQUAL(2, bins[2]);

    free(bins);
    bins = NULL;

    /* 
     * > 2.0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     *  3	7	12.0	13.0	0	0
     *  3	8	13.0	14.0	0	0
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     8,
                                                     0.0,
                                                     2.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(7,num_query_bins);
    TEST_ASSERT_EQUAL(3, bins[0]);
    TEST_ASSERT_EQUAL(4, bins[1]);
    TEST_ASSERT_EQUAL(5, bins[2]);
    TEST_ASSERT_EQUAL(6, bins[3]);
    TEST_ASSERT_EQUAL(7, bins[4]);
    TEST_ASSERT_EQUAL(8, bins[5]);
    TEST_ASSERT_EQUAL(9, bins[6]);

    free(bins);
    bins = NULL;

    /*
     * >= 2
     *  3	2	2.0	3.0	0	0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     *  3	7	12.0	13.0	0	0
     *  3	8	13.0	14.0	0	0
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     16,
                                                     0.0,
                                                     2.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(8,num_query_bins);
    TEST_ASSERT_EQUAL(2, bins[0]);
    TEST_ASSERT_EQUAL(3, bins[1]);
    TEST_ASSERT_EQUAL(4, bins[2]);
    TEST_ASSERT_EQUAL(5, bins[3]);
    TEST_ASSERT_EQUAL(6, bins[4]);
    TEST_ASSERT_EQUAL(7, bins[5]);
    TEST_ASSERT_EQUAL(8, bins[6]);
    TEST_ASSERT_EQUAL(9, bins[7]);

    free(bins);
    bins = NULL;

    /* 
     * > when hiting lo bi
     *  3	0	0.0	1.0	1	0
     *  3	1	1.0	2.0	0	0
     *  3	2	2.0	3.0	0	0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     *  3	7	12.0	13.0	0	0
     *  3	8	13.0	14.0	0	0
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     8,
                                                     0.0,
                                                     -2.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(10,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);
    TEST_ASSERT_EQUAL(1, bins[1]);
    TEST_ASSERT_EQUAL(2, bins[2]);
    TEST_ASSERT_EQUAL(3, bins[3]);
    TEST_ASSERT_EQUAL(4, bins[4]);
    TEST_ASSERT_EQUAL(5, bins[5]);
    TEST_ASSERT_EQUAL(6, bins[6]);
    TEST_ASSERT_EQUAL(7, bins[7]);
    TEST_ASSERT_EQUAL(8, bins[8]);
    TEST_ASSERT_EQUAL(9, bins[9]);

    free(bins);
    bins = NULL;

    /*
     * > when hiting hi bi
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     8,
                                                     2.0,
                                                     20.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(9, bins[0]);

    free(bins);
    bins = NULL;

    /*
     * < 10 
     * > 2
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     2+8,
                                                     10.0, // <
                                                     2.0, // >
                                                     &bins);
    TEST_ASSERT_EQUAL(4,num_query_bins);
    TEST_ASSERT_EQUAL(3, bins[0]);
    TEST_ASSERT_EQUAL(4, bins[1]);
    TEST_ASSERT_EQUAL(5, bins[2]);
    TEST_ASSERT_EQUAL(6, bins[3]);


    free(bins);
    bins = NULL;

    /* 
     * <= 10 
     * > 2
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     4+8,
                                                     10.0, // <
                                                     2.0, // >
                                                     &bins);
    TEST_ASSERT_EQUAL(4,num_query_bins);
    TEST_ASSERT_EQUAL(3, bins[0]);
    TEST_ASSERT_EQUAL(4, bins[1]);
    TEST_ASSERT_EQUAL(5, bins[2]);
    TEST_ASSERT_EQUAL(6, bins[3]);

    free(bins);
    bins = NULL;

    /* 
     * <= 10
     * >= 2
     *  3	2	2.0	3.0	0	0
     *  3	3	3.0	4.0	0	0
     *  3	4	4.0	6.0	0	0
     *  3	5	6.0	8.0	0	0
     *  3	6	8.0	12.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     4+16,
                                                     10.0, // <
                                                     2.0, // >
                                                     &bins);
    TEST_ASSERT_EQUAL(5,num_query_bins);
    TEST_ASSERT_EQUAL(2, bins[0]);
    TEST_ASSERT_EQUAL(3, bins[1]);
    TEST_ASSERT_EQUAL(4, bins[2]);
    TEST_ASSERT_EQUAL(5, bins[3]);
    TEST_ASSERT_EQUAL(6, bins[4]);

    free(bins);
    bins = NULL;

    /*  
     * Range when hitting lo bin on both
     * <= -1
     * >= -2
     *  3	0	0.0	1.0	1	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     4+16,
                                                     -1.0, // <
                                                     -2.0, // >
                                                     &bins);
    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0, bins[0]);

    free(bins);
    bins = NULL;

    /*  
     * Range when hitting hi bin on both
     * <= 100 
     * >= 50 
     *  3	9	14.0	0.0	0	1
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "AC",
                                                     ac_rowid,
                                                     num_bins,
                                                     ac_md_type,
                                                     4+16,
                                                     100.0, // <
                                                     50.0, // >
                                                     &bins);
    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(9, bins[0]);

    free(bins);
    bins = NULL;
    //}}}

    int lof_md_type;
    int lof_rowid = get_variant_metadata_bin_info(variant_db_name,
                                                  "LOF",
                                                  &offset,
                                                  &num_bins,
                                                  &lof_md_type,
                                                  &source_file);

    //{{{ TEST LOF
    /*
     * select * from variant_md_bins where md_rowid==1;
     * 1	0	0.0	1.0	0	0
     * 1	1	1.0	2.0	0	0
     */

    // == 5
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     1,
                                                     5.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;

    /* 
     * == 0
     * 1	0	0.0	1.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     1,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0,bins[0]);

    free(bins);
    bins = NULL;


    /* 
     * == 1
     * 1	1	1.0	2.0	0	0
     */
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     1,
                                                     1.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(1,bins[0]);

    free(bins);
    bins = NULL;


    // < 0
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     2,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;

    // > 1
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     8,
                                                     0.0,
                                                     1.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;


   
     /* 
      * <= 0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     4,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;


     /* 
      * <= 1
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     4,
                                                     1.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;


     /* 
      * < 2
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     2,
                                                     2.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;


    // > 1
    num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     8,
                                                     0.0,
                                                     1.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;

   
    /* 
     * >= 1
     */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     16,
                                                     0.0,
                                                     1.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;



     /* 
      * >= 0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     16,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;


     /* 
      * > 0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     8,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;

     /* 
      * > -1 
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     8,
                                                     0.0,
                                                     -1.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;

 
     /* 
      * >= -1 
      * 1	0	0.0	1.0	0	0
      * 1	1	1.0	2.0	0	0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     16,
                                                     0.0,
                                                     -1.0,
                                                     &bins);
    TEST_ASSERT_EQUAL(0,num_query_bins);
    free(bins);
    bins = NULL;

     /* 
      * != 1 
      * 1	0	0.0	1.0	0	0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     32,
                                                     1.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(0,bins[0]);

    free(bins);
    bins = NULL;

     /* 
      * != 0 
      * 1	1	1.0	2.0	0	0
      */
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     32,
                                                     0.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(1,num_query_bins);
    TEST_ASSERT_EQUAL(1,bins[0]);

    free(bins);
    bins = NULL;


     // > 0 < 1
     num_query_bins = get_variant_metadata_query_bins(variant_db_name,
                                                     "LOF",
                                                     lof_rowid,
                                                     num_bins,
                                                     lof_md_type,
                                                     3+8,
                                                     1.0,
                                                     0.0,
                                                     &bins);

    TEST_ASSERT_EQUAL(0,num_query_bins);

    free(bins);
    bins = NULL;
    //}}}

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);
}
//}}}

void test_pares_q(void)
{
    struct gqt_query q;
    
    int r = parse_q("count(HET)", &q);
}

void test_get_variant_metadata_wah_bitmap_in_place(void)
{
    char *gt_file_name = ".test_get_query_bins.gqt",
         *bim_file_name = ".test_get_query_bins.bim",
         *vid_file_name = ".test_get_query_bins.vid",
         *tmp_dir = ".",
         *variant_index_file_name = ".test_get_query_bins.gqtv",
         *variant_db_name = ".test_get_query_bins.gqtv.db";

    char *full_cmd = "test_get_variant_metadata_wah_bitmap_in_place";

    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);

    //{{{ Index LOF, AF, and AC
    uint32_t r = convert_file_by_name_bcf_to_wahbm_bim(BCF_FILE,
                                                       NUM_INDS,
                                                       NUM_VARS,
                                                       gt_file_name,
                                                       bim_file_name,
                                                       vid_file_name,
                                                       tmp_dir,
                                                       full_cmd);


    char *LOF_field = "LOF";
    uint32_t num_to_test = 20, lof_num_bins = 50;
    void *bin_range_lo = NULL, *bin_range_hi = NULL;
    int less_than_bin = -1, greater_than_bin = -1;

    int type = index_variant_metadata(BCF_FILE,
                                      vid_file_name,
                                      variant_db_name,
                                      variant_index_file_name,
                                      NUM_VARS,
                                      LOF_field,
                                      num_to_test,
                                      &lof_num_bins,
                                      &bin_range_lo,
                                      &bin_range_hi,
                                      &less_than_bin,
                                      &greater_than_bin);

    char *AF_field = "AF";
    num_to_test = 20;
    uint32_t af_num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AF_field,
                                  num_to_test,
                                  &af_num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    char *AC_field = "AC";
    num_to_test = 20;
    uint32_t ac_num_bins = 50;
    free(bin_range_lo);
    free(bin_range_hi);
    bin_range_lo = NULL;
    bin_range_hi = NULL;
    less_than_bin = -1;
    greater_than_bin = -1;

    type = index_variant_metadata(BCF_FILE,
                                  vid_file_name,
                                  variant_db_name,
                                  variant_index_file_name,
                                  NUM_VARS,
                                  AC_field,
                                  num_to_test,
                                  &ac_num_bins,
                                  &bin_range_lo,
                                  &bin_range_hi,
                                  &less_than_bin,
                                  &greater_than_bin);

    char *source_file;
    uint64_t offset;
    uint32_t num_bins;
    int af_md_type;
    int af_rowid = get_variant_metadata_bin_info(variant_db_name,
                                                 "AF",
                                                 &offset,
                                                 &num_bins,
                                                 &af_md_type,
                                                 &source_file);
    //}}}

    int lof_md_type;
    int lof_rowid = get_variant_metadata_bin_info(variant_db_name,
                                                  "LOF",
                                                  &offset,
                                                  &num_bins,
                                                  &lof_md_type,
                                                  &source_file);

    //fprintf(stderr, "%llu\n", offset);

    FILE *f = fopen(variant_index_file_name, "rb");
    fseek(f, offset, SEEK_SET);
    uint32_t *bin_sizes = (uint32_t *)malloc(num_bins*sizeof(uint32_t));
    r = fread(bin_sizes, sizeof(uint32_t), num_bins, f); 

    uint32_t i;
    /*
    for (i = 0; i < num_bins; ++i)
        fprintf(stderr, "%u\n", bin_sizes[i]);
    */



    remove(gt_file_name);
    remove(bim_file_name);
    remove(vid_file_name);
    remove(variant_index_file_name);
    remove(variant_db_name);
}

//{{{void test_put_it_all_together(void)
void test_put_it_all_together(void)
{
    struct bcf_file bcf_f = init_bcf_file(BCF_FILE);
    float *LOF;
    float LOF_lo[2] = {0,1};
    float LOF_hi[2] = {1,2};
    float missing = 0;
    int type = 0;

    int r = get_variant_metadata(&bcf_f,
                                 NUM_VARS,
                                 "LOF",
                                 (void **)&LOF,
                                 &missing,
                                 &type);

    int to_test = 20, num_bins = 2;

    float *bin_range_lo = NULL, *bin_range_hi = NULL;

    uint32_t actual_num_bins = float_equal_freq_binning(LOF,
                                                      to_test,
                                                      num_bins,
                                                      &bin_range_lo,
                                                      &bin_range_hi);
    TEST_ASSERT_EQUAL(num_bins, actual_num_bins);

    uint32_t i,j;
    for (i = 0; i < actual_num_bins; ++i) {
        TEST_ASSERT_EQUAL(LOF_lo[i], bin_range_lo[i]);
        TEST_ASSERT_EQUAL(LOF_hi[i], bin_range_hi[i]);
    }

    close_bcf_file(&bcf_f);

    uint32_t **ws = NULL;
    uint32_t *w_lens = NULL;

    uint32_t bit_array_lens = float_to_wah_bitmap(LOF,
                                                NUM_VARS,
                                                bin_range_lo,
                                                bin_range_hi,
                                                2,
                                                0,
                                                0,
                                                2,
                                                &ws,
                                                &w_lens);
}
///}}}

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
