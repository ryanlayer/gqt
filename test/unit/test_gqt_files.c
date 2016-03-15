#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <htslib/knetfile.h>

#include "unity.h"
#include "wahbm.h"
#include "parse_q.h"
#include "genotq.h"
#include "vid.h"
#include "bim.h"
#include "ubin.h"

char *BCF_FILE = "../data/diff_gts.bcf";
uint32_t NUM_INDS = 10;
uint32_t NUM_VARS = 43;

void setUp(void) { }
void tearDown(void) { }


//{{{void test_get_wahbm_bitmap(void)
void test_get_wahbm_bitmap(void)
{
    char *file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 60;
    uint32_t num_samples = 3;

    struct wahbm_file *wf = new_wahbm_file(file_name,
                                           full_cmd,
                                           num_variants,
                                           num_samples);



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

    append_ubin_to_wahbm_file(wf, 0, ubin_0, ubin_len_0);

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, 60, &ubin_1);

    append_ubin_to_wahbm_file(wf, 1, ubin_1, ubin_len_1);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, 60, &ubin_2);

    append_ubin_to_wahbm_file(wf, 2, ubin_2, ubin_len_2);

    write_record_offsets(wf);

    destroy_wahbm_file(wf);

    struct wahbm_file *o = open_wahbm_file(file_name);

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

//{{{ void test_open_wahbm_file(void)
void test_open_wahbm_file(void)
{
    char *file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 60;
    uint32_t num_samples = 3;

    struct wahbm_file *w = new_wahbm_file(file_name,
                                          full_cmd,
                                          num_variants,
                                          num_samples);



    char *plt_0 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_1 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_2 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";


    uint32_t *ubin_0;
    uint32_t ubin_len_0 = plt_line_to_packed_ints(plt_0, 60, &ubin_0);

    uint32_t *wah_0;
    uint32_t *wah_sizes_0;
    uint32_t wah_len_0 = ubin_to_bitmap_wah(ubin_0,
                                            ubin_len_0,
                                            60,
                                            &wah_0,
                                            &wah_sizes_0);

    append_wahbm_to_wahbm_file(w, 0, wah_0, wah_len_0, wah_sizes_0);

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, 60, &ubin_1);

    uint32_t *wah_1;
    uint32_t *wah_sizes_1;
    uint32_t wah_len_1 = ubin_to_bitmap_wah(ubin_1,
                                            ubin_len_1,
                                            60,
                                            &wah_1,
                                            &wah_sizes_1);

    append_wahbm_to_wahbm_file(w, 1, wah_1, wah_len_1, wah_sizes_1);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, 60, &ubin_2);

    uint32_t *wah_2;
    uint32_t *wah_sizes_2;
    uint32_t wah_len_2 = ubin_to_bitmap_wah(ubin_2,
                                            ubin_len_2,
                                            60,
                                            &wah_2,
                                            &wah_sizes_2);

    append_wahbm_to_wahbm_file(w, 2, wah_2, wah_len_2, wah_sizes_2);

    write_record_offsets(w);

    destroy_wahbm_file(w);

    struct wahbm_file *o = open_wahbm_file(file_name);

    TEST_ASSERT_EQUAL('G', o->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', o->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', o->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('w', o->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), o->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), o->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), o->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), o->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, o->gqt_header->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), o->gqt_header->id_hash);
    TEST_ASSERT_EQUAL(num_variants, o->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(num_samples, o->gqt_header->num_samples);

    TEST_ASSERT_EQUAL(sizeof(struct gqt_file_header) + 
                            4*num_samples*sizeof(uint64_t), 
                      o->header_offset);



    /*
    FILE *f = fopen(file_name, "rb");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", file_name);

    if (fseek(f, sizeof(struct gqt_file_header), SEEK_SET) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            file_name);

    uint64_t *R = (uint64_t *) malloc (4*num_samples*sizeof(uint64_t));

    if (fread(R, sizeof(uint64_t), 4*num_samples, f) != 4*num_samples)
        err(EX_IOERR, "Read error '%s'", file_name);

    uint32_t i;
    for (i = 0; i < num_samples*4; ++i) 
        TEST_ASSERT_EQUAL(running_total[i], R[i]);


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

    uint32_t V;

    for (i = 0; i < num_samples*4; ++i) {
        if (fread(&V, sizeof(uint32_t), 1, f) != 1)
            err(EX_IOERR, "Read error '%s'", file_name);
        TEST_ASSERT_EQUAL( A[i], V);
    }

    fclose(f);
    */
}
//}}}

//{{{void test_append_ubin_to_wahbm_file(void)
void test_append_ubin_to_wahbm_file(void)
{

    char *file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 60;
    uint32_t num_samples = 3;

    struct wahbm_file *wf = new_wahbm_file(file_name,
                                           full_cmd,
                                           num_variants,
                                           num_samples);



    char *plt_0 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_1 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_2 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";


    uint32_t *ubin_0;
    uint32_t ubin_len_0 = plt_line_to_packed_ints(plt_0, 60, &ubin_0);

    append_ubin_to_wahbm_file(wf, 0, ubin_0, ubin_len_0);

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, 60, &ubin_1);

    append_ubin_to_wahbm_file(wf, 1, ubin_1, ubin_len_1);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, 60, &ubin_2);

    append_ubin_to_wahbm_file(wf, 2, ubin_2, ubin_len_2);

    write_record_offsets(wf);

    destroy_wahbm_file(wf);

    FILE *f = fopen(file_name, "rb");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", file_name);

    if (fseek(f, sizeof(struct gqt_file_header), SEEK_SET) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            file_name);

    uint64_t *R = (uint64_t *) malloc (4*num_samples*sizeof(uint64_t));

    if (fread(R, sizeof(uint64_t), 4*num_samples, f) != 4*num_samples)
        err(EX_IOERR, "Read error '%s'", file_name);

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

    uint32_t i, V;

    for (i = 0; i < num_samples*4; ++i) {
        if (fread(&V, sizeof(uint32_t), 1, f) != 1)
            err(EX_IOERR, "Read error '%s'", file_name);
        TEST_ASSERT_EQUAL( A[i], V);
    }

    fclose(f);
}
//}}}

//{{{void test_append_wahbm_to_wahbm_file(void)
void test_append_wahbm_to_wahbm_file(void)
{

    char *file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 60;
    uint32_t num_samples = 3;

    struct wahbm_file *w = new_wahbm_file(file_name,
                                          full_cmd,
                                          num_variants,
                                          num_samples);



    char *plt_0 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";

    char *plt_1 = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";

    char *plt_2 = "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 "
                  "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0";


    uint32_t *ubin_0;
    uint32_t ubin_len_0 = plt_line_to_packed_ints(plt_0, 60, &ubin_0);

    uint32_t *wah_0;
    uint32_t *wah_sizes_0;
    uint32_t wah_len_0 = ubin_to_bitmap_wah(ubin_0,
                                            ubin_len_0,
                                            60,
                                            &wah_0,
                                            &wah_sizes_0);

    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111111"),
                       wah_0[0]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111100"),
                       wah_0[1]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_0[2]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_0[3]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_0[4]);



    append_wahbm_to_wahbm_file(w, 0, wah_0, wah_len_0, wah_sizes_0);

    uint64_t *running_total = (uint64_t *)
            malloc(4*num_samples*sizeof(uint64_t));
    running_total[0] = wah_sizes_0[0];
    TEST_ASSERT_EQUAL(running_total[0], w->record_offsets[0]);

    running_total[1] = running_total[0] + wah_sizes_0[1];
    TEST_ASSERT_EQUAL(running_total[1], w->record_offsets[1]);

    running_total[2] = running_total[1] + wah_sizes_0[2];
    TEST_ASSERT_EQUAL(running_total[2], w->record_offsets[2]);

    running_total[3] = running_total[2] + wah_sizes_0[3];
    TEST_ASSERT_EQUAL(running_total[3], w->record_offsets[3]);

    uint32_t *ubin_1;
    uint32_t ubin_len_1 = plt_line_to_packed_ints(plt_1, 60, &ubin_1);

    uint32_t *wah_1;
    uint32_t *wah_sizes_1;
    uint32_t wah_len_1 = ubin_to_bitmap_wah(ubin_1,
                                            ubin_len_1,
                                            60,
                                            &wah_1,
                                            &wah_sizes_1);

    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[0]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111111"),
                       wah_1[1]);
    TEST_ASSERT_EQUAL( bin_char_to_int("01111111111111111111111111111100"),
                       wah_1[2]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[3]);
    TEST_ASSERT_EQUAL( bin_char_to_int("10000000000000000000000000000010"),
                       wah_1[4]);

    append_wahbm_to_wahbm_file(w, 1, wah_1, wah_len_1, wah_sizes_1);

    ///// MOVE RUNNING TOTAL TO ARRAY!
    running_total[4] = running_total[3] + wah_sizes_1[0];
    TEST_ASSERT_EQUAL(running_total[4], w->record_offsets[4]);

    running_total[5] = running_total[4] + wah_sizes_1[1];
    TEST_ASSERT_EQUAL(running_total[5], w->record_offsets[5]);

    running_total[6] = running_total[5] + wah_sizes_1[2];
    TEST_ASSERT_EQUAL(running_total[6], w->record_offsets[6]);

    running_total[7] = running_total[6] + wah_sizes_1[3];
    TEST_ASSERT_EQUAL(running_total[7], w->record_offsets[7]);

    uint32_t *ubin_2;
    uint32_t ubin_len_2 = plt_line_to_packed_ints(plt_2, 60, &ubin_2);

    uint32_t *wah_2;
    uint32_t *wah_sizes_2;
    uint32_t wah_len_2 = ubin_to_bitmap_wah(ubin_2,
                                            ubin_len_2,
                                            60,
                                            &wah_2,
                                            &wah_sizes_2);

    append_wahbm_to_wahbm_file(w, 2, wah_2, wah_len_2, wah_sizes_2);

    running_total[8] = running_total[7] + wah_sizes_2[0];
    TEST_ASSERT_EQUAL(running_total[8], w->record_offsets[8]);

    running_total[9] = running_total[8] + wah_sizes_2[1];
    TEST_ASSERT_EQUAL(running_total[9], w->record_offsets[9]);

    running_total[10] = running_total[9] + wah_sizes_2[2];
    TEST_ASSERT_EQUAL(running_total[10], w->record_offsets[10]);

    running_total[11] = running_total[10] + wah_sizes_2[3];
    TEST_ASSERT_EQUAL(running_total[11], w->record_offsets[11]);

    write_record_offsets(w);

    destroy_wahbm_file(w);

    FILE *f = fopen(file_name, "rb");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", file_name);

    if (fseek(f, sizeof(struct gqt_file_header), SEEK_SET) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            file_name);

    uint64_t *R = (uint64_t *) malloc (4*num_samples*sizeof(uint64_t));

    if (fread(R, sizeof(uint64_t), 4*num_samples, f) != 4*num_samples)
        err(EX_IOERR, "Read error '%s'", file_name);

    uint32_t i;
    for (i = 0; i < num_samples*4; ++i) 
        TEST_ASSERT_EQUAL(running_total[i], R[i]);


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

    uint32_t V;

    for (i = 0; i < num_samples*4; ++i) {
        if (fread(&V, sizeof(uint32_t), 1, f) != 1)
            err(EX_IOERR, "Read error '%s'", file_name);
        TEST_ASSERT_EQUAL( A[i], V);
    }

    fclose(f);
}
//}}}

//{{{void test_new_gqt_file(void)
void test_new_gqt_file(void)
{
    char *file_name = "test_gqt";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 4;
    uint32_t num_samples = 4;

    struct wahbm_file *w = new_wahbm_file(file_name,
                                          full_cmd,
                                          num_variants,
                                          num_samples);

    destroy_wahbm_file(w);

    FILE *f = fopen(file_name, "rb");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", file_name);

    size_t r;

    char c;

    //char marker[3]; // "GQT"
    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('G', c);

    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('Q', c);

    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('T', c);

    //char type; // g gqt, v vid, b bim
    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('w', c);


    uint32_t u32;
    //uint32_t major, minor, revision, build;
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), u32);

    //uint32_t magic; //0x11223344
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(0x11223344, u32);

    unsigned long h;
    //unsigned long id_hash;// used to varify the files were created together.
    if (fread(&h, sizeof(unsigned long), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), h);

  
    //uint32_t num_variants, num_samples;
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(num_variants, u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(num_samples, u32);

    uint32_t m[20];
    //uint32_t more[20];
    if (fread(&m, sizeof(uint32_t), 20, f) != 20)
        err(EX_IOERR, "Read error '%s'", file_name);
    for (u32 = 0; u32 < 20; u32++)
        TEST_ASSERT_EQUAL(0, m[u32]);

    uint32_t u64;
    //uint64_t *record_offsets;
    for (u32 = 0; u32 < 4*num_samples; u32++) {
        if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
            err(EX_IOERR, "Read error '%s'", file_name);
        TEST_ASSERT_EQUAL(0, u64);
    }

    // we should be at the end of the file and this read should fail
    TEST_ASSERT_EQUAL(0, fread(&u64, sizeof(uint64_t), 1, f));
    TEST_ASSERT_EQUAL(1, feof(f));

    fclose(f);

    remove(file_name);
}
//}}}

//{{{ void test_seek_bim_to_data(void)
void test_seek_bim_to_data(void)
{
    char *file_name = "test_bim";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 4;
    uint32_t num_samples = 4;
    uint64_t u_size = 1;
    uint64_t c_size = 2;
    uint64_t h_size = 3;
    uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
    md_line_lens[0] = 1;
    md_line_lens[1] = 2;
    md_line_lens[2] = 3;
    md_line_lens[3] = 4;

    struct bim_file *b = new_bim_file(file_name,
                                      full_cmd,
                                      num_variants,
                                      num_samples,
                                      u_size,
                                      c_size,
                                      h_size,
                                      md_line_lens);

    uint32_t test_data[4] = {10, 11, 12, 13};
    if (fwrite(test_data,
               sizeof(uint64_t),
               4,
               b->file.local) != 4)
        err(EX_IOERR, "Error writing to BIM file \"%s\"", file_name);

    destroy_bim_file(b);

    struct bim_file *b1 = open_bim_file(file_name);

    TEST_ASSERT_EQUAL('G', b1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', b1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', b1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('b', b1->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), b1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), b1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), b1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), b1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, b1->gqt_header->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), b1->gqt_header->id_hash);
    TEST_ASSERT_EQUAL(num_variants, b1->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(num_samples, b1->gqt_header->num_samples);

    uint32_t i;
    for (i = 0; i < 20; i++)
        TEST_ASSERT_EQUAL(0, b1->gqt_header->more[i]);

    TEST_ASSERT_EQUAL(u_size, b1->bim_header->u_size);
    TEST_ASSERT_EQUAL(c_size, b1->bim_header->c_size);
    TEST_ASSERT_EQUAL(h_size, b1->bim_header->h_size);

    TEST_ASSERT_EQUAL(1, b1->bim_header->md_line_lens[0]);
    TEST_ASSERT_EQUAL(2, b1->bim_header->md_line_lens[1]);
    TEST_ASSERT_EQUAL(3, b1->bim_header->md_line_lens[2]);
    TEST_ASSERT_EQUAL(4, b1->bim_header->md_line_lens[3]);

    seek_bim_to_data(b1);

    uint32_t read_data_0[4];
    //if (fread(read_data_0, sizeof(uint32_t), 4, b1->file.remote) != 4)
    if (knet_read(b1->file.remote, 
                  read_data_0, 4 * sizeof(uint32_t)) != 4 * sizeof(uint32_t))
        err(EX_IOERR, "Read error '%s'", b1->file_name);
    TEST_ASSERT_EQUAL(10, read_data_0[0]);
    TEST_ASSERT_EQUAL(11, read_data_0[1]);
    TEST_ASSERT_EQUAL(12, read_data_0[2]);
    TEST_ASSERT_EQUAL(13, read_data_0[3]);

    seek_bim_to_data(b1);

    uint32_t read_data_1[4];
    //if (fread(read_data_1, sizeof(uint32_t), 4, b1->file.remote) != 4)
    if (knet_read(b1->file.remote,
                  read_data_1, 
                  4 * sizeof(uint32_t)) != 4 * sizeof(uint32_t))
        err(EX_IOERR, "Read error '%s'", b1->file_name);
    TEST_ASSERT_EQUAL(10, read_data_1[0]);
    TEST_ASSERT_EQUAL(11, read_data_1[1]);
    TEST_ASSERT_EQUAL(12, read_data_1[2]);
    TEST_ASSERT_EQUAL(13, read_data_1[3]);

    destroy_bim_file(b1);

    remove(file_name);

}
//}}}

//{{{ void test_upate_bim_header(void)
void test_upate_bim_header(void)
{
    char *file_name = "test_bim";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 4;
    uint32_t num_samples = 4;
    uint64_t u_size = 1;
    uint64_t c_size = 2;
    uint64_t h_size = 3;
    uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
    md_line_lens[0] = 1;
    md_line_lens[1] = 2;
    md_line_lens[2] = 3;
    md_line_lens[3] = 4;

    struct bim_file *b = new_bim_file(file_name,
                                      full_cmd,
                                      num_variants,
                                      num_samples,
                                      u_size,
                                      c_size,
                                      h_size,
                                      md_line_lens);

    //destroy_bim_file(b);

    //struct bim_file *b1 = open_bim_file(file_name);
    struct bim_file *b1 = b;

    TEST_ASSERT_EQUAL('G', b1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', b1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', b1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('b', b1->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), b1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), b1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), b1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), b1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, b1->gqt_header->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), b1->gqt_header->id_hash);
    TEST_ASSERT_EQUAL(num_variants, b1->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(num_samples, b1->gqt_header->num_samples);

    uint32_t i;
    for (i = 0; i < 20; i++)
        TEST_ASSERT_EQUAL(0, b1->gqt_header->more[i]);

    TEST_ASSERT_EQUAL(u_size, b1->bim_header->u_size);
    TEST_ASSERT_EQUAL(c_size, b1->bim_header->c_size);
    TEST_ASSERT_EQUAL(h_size, b1->bim_header->h_size);

    TEST_ASSERT_EQUAL(1, b1->bim_header->md_line_lens[0]);
    TEST_ASSERT_EQUAL(2, b1->bim_header->md_line_lens[1]);
    TEST_ASSERT_EQUAL(3, b1->bim_header->md_line_lens[2]);
    TEST_ASSERT_EQUAL(4, b1->bim_header->md_line_lens[3]);

    update_bim_file_header(10, 11, 12, b1);

    TEST_ASSERT_EQUAL(10, b1->bim_header->u_size);
    TEST_ASSERT_EQUAL(11, b1->bim_header->c_size);
    TEST_ASSERT_EQUAL(12, b1->bim_header->h_size);

    //destroy_bim_file(b1);

    struct bim_file *b2 = b;

    TEST_ASSERT_EQUAL('G', b2->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', b2->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', b2->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('b', b2->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), b2->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), b2->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), b2->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), b2->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, b2->gqt_header->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), b2->gqt_header->id_hash);
    TEST_ASSERT_EQUAL(num_variants, b2->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(num_samples, b2->gqt_header->num_samples);

    for (i = 0; i < 20; i++)
        TEST_ASSERT_EQUAL(0, b2->gqt_header->more[i]);

    TEST_ASSERT_EQUAL(10, b2->bim_header->u_size);
    TEST_ASSERT_EQUAL(11, b2->bim_header->c_size);
    TEST_ASSERT_EQUAL(12, b2->bim_header->h_size);

    TEST_ASSERT_EQUAL(1, b2->bim_header->md_line_lens[0]);
    TEST_ASSERT_EQUAL(2, b2->bim_header->md_line_lens[1]);
    TEST_ASSERT_EQUAL(3, b2->bim_header->md_line_lens[2]);
    TEST_ASSERT_EQUAL(4, b2->bim_header->md_line_lens[3]);

    destroy_bim_file(b2);

    remove(file_name);
}
//}}}

//{{{void test_open_bim_file(void)
void test_open_bim_file(void)
{
    char *file_name = "test_bim";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 4;
    uint32_t num_samples = 4;
    uint64_t u_size = 1;
    uint64_t c_size = 2;
    uint64_t h_size = 3;
    uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
    md_line_lens[0] = 1;
    md_line_lens[1] = 2;
    md_line_lens[2] = 3;
    md_line_lens[3] = 4;

    struct bim_file *b = new_bim_file(file_name,
                                      full_cmd,
                                      num_variants,
                                      num_samples,
                                      u_size,
                                      c_size,
                                      h_size,
                                      md_line_lens);

    uint64_t header_size = sizeof(struct gqt_file_header) +
                           sizeof(uint64_t) * 3 +
                           sizeof(uint64_t) * 4;

    TEST_ASSERT_EQUAL(header_size,b->data_start);

    destroy_bim_file(b);

    struct bim_file *b1 = open_bim_file(file_name);

    TEST_ASSERT_EQUAL('G', b1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', b1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', b1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('b', b1->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), b1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), b1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), b1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), b1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, b1->gqt_header->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), b1->gqt_header->id_hash);
    TEST_ASSERT_EQUAL(num_variants, b1->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(num_samples, b1->gqt_header->num_samples);
    TEST_ASSERT_EQUAL(header_size,b1->data_start);

    uint32_t i;
    for (i = 0; i < 20; i++)
        TEST_ASSERT_EQUAL(0, b1->gqt_header->more[i]);

    TEST_ASSERT_EQUAL(u_size, b1->bim_header->u_size);
    TEST_ASSERT_EQUAL(c_size, b1->bim_header->c_size);
    TEST_ASSERT_EQUAL(h_size, b1->bim_header->h_size);

    TEST_ASSERT_EQUAL(1, b1->bim_header->md_line_lens[0]);
    TEST_ASSERT_EQUAL(2, b1->bim_header->md_line_lens[1]);
    TEST_ASSERT_EQUAL(3, b1->bim_header->md_line_lens[2]);
    TEST_ASSERT_EQUAL(4, b1->bim_header->md_line_lens[3]);

    destroy_bim_file(b1);

    remove(file_name);
}
//}}}

//{{{void test_new_bim_file(void)
void test_new_bim_file(void)
{
    char *file_name = "test_bim";
    char *full_cmd = "gqt convert bcf -i bcf";
    uint32_t num_variants = 4;
    uint32_t num_samples = 4;
    uint64_t u_size = 1;
    uint64_t c_size = 2;
    uint64_t h_size = 3;
    uint64_t *md_line_lens = (uint64_t *) malloc(4*sizeof(uint64_t));
    md_line_lens[0] = 1;
    md_line_lens[1] = 2;
    md_line_lens[2] = 3;
    md_line_lens[3] = 4;

    struct bim_file *b = new_bim_file(file_name,
                                      full_cmd,
                                      num_variants,
                                      num_samples,
                                      u_size,
                                      c_size,
                                      h_size,
                                      md_line_lens);

    destroy_bim_file(b);

    FILE *f = fopen(file_name, "rb");
    if (!f)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", file_name);

    size_t r;

    char c;

    //char marker[3]; // "GQT"
    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('G', c);

    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('Q', c);

    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('T', c);

    //char type; // g gqt, v vid, b bim
    if (fread(&c, sizeof(char), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL('b', c);


    uint32_t u32;
    //uint32_t major, minor, revision, build;
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), u32);

    //uint32_t magic; //0x11223344
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(0x11223344, u32);

    unsigned long h;
    //unsigned long id_hash;// used to varify the files were created together.
    if (fread(&h, sizeof(unsigned long), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), h);

  
    //uint32_t num_variants, num_samples;
    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(num_variants, u32);

    if (fread(&u32, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(num_samples, u32);

    uint32_t m[20];
    //uint32_t more[20];
    if (fread(&m, sizeof(uint32_t), 20, f) != 20)
        err(EX_IOERR, "Read error '%s'", file_name);
    for (u32 = 0; u32 < 20; u32++)
        TEST_ASSERT_EQUAL(0, m[u32]);

    uint32_t u64;
    //uint64_t u_size, c_size, h_size; 
    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(u_size, u64);

    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(c_size, u64);

    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(h_size, u64);

    //uint64_t *md_line_lens;
    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(1, u64);

    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(2, u64);

    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(3, u64);

    if (fread(&u64, sizeof(uint64_t), 1, f) != 1)
        err(EX_IOERR, "Read error '%s'", file_name);
    TEST_ASSERT_EQUAL(4, u64);

    fclose(f);

    remove(file_name);
}
//}}}

//{{{void test_new_bim_file_header(void)
void test_new_bim_file_header(void)
{
    uint64_t u_size = 1;
    uint64_t c_size = 2;
    uint64_t h_size = 3;
    uint64_t md_line_lens[4] = {1, 2, 3, 4};

    struct bim_file_header *h = new_bim_file_header(u_size,
                                                    c_size,
                                                    h_size,
                                                    md_line_lens);
    TEST_ASSERT_EQUAL(u_size, h->u_size);
    TEST_ASSERT_EQUAL(c_size, h->c_size);
    TEST_ASSERT_EQUAL(h_size, h->h_size);

    TEST_ASSERT_EQUAL(md_line_lens[0], h->md_line_lens[0]);
    TEST_ASSERT_EQUAL(md_line_lens[1], h->md_line_lens[1]);
    TEST_ASSERT_EQUAL(md_line_lens[2], h->md_line_lens[2]);
    TEST_ASSERT_EQUAL(md_line_lens[3], h->md_line_lens[3]);
}
//}}}

//{{{void test_new_gqt_file_header(void)
void test_new_gqt_file_header(void)
{
    char *full_cmd = "gqt convert bcf -i bcf";
    struct gqt_file_header *h = new_gqt_file_header('v',
                                                    full_cmd,
                                                    NUM_VARS,
                                                    NUM_INDS);

    // cannot start with a number
    TEST_ASSERT_EQUAL('v', h->type);

    TEST_ASSERT_EQUAL('G', h->marker[0]);
    TEST_ASSERT_EQUAL('Q', h->marker[1]);
    TEST_ASSERT_EQUAL('T', h->marker[2]);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), h->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), h->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), h->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), h->build);
    TEST_ASSERT_EQUAL(0x11223344, h->magic);
    TEST_ASSERT_EQUAL(hash_cmd(full_cmd), h->id_hash);

    uint32_t i;
    for (i = 0; i < 20; ++i)
        TEST_ASSERT_EQUAL(0, h->more[i]);
}
//}}}

//{{{void test_new_and_open_vid_file(void)
void test_new_and_open_vid_file(void)
{
    char *full_cmd = "gqt convert bcf -i bcf";
    struct vid_file *v0 = new_vid_file("test_vid",
                                       full_cmd,
                                       NUM_VARS,
                                       NUM_INDS);

    TEST_ASSERT_EQUAL('G', v0->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', v0->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', v0->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), v0->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), v0->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), v0->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), v0->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, v0->gqt_header->magic);

    destroy_vid_file(v0);

    struct vid_file *v1 = open_vid_file("test_vid");

    TEST_ASSERT_EQUAL('G', v1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', v1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', v1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), v1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), v1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), v1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), v1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, v1->gqt_header->magic);

    destroy_vid_file(v1);

    remove("test_vid");
}
//}}}

//{{{void test_write_lod_vid_file(void)
void test_write_lod_vid_file(void)
{
    char *full_cmd = "gqt convert bcf -i test.bcf";
    struct vid_file *v0 = new_vid_file("test_vid",
                                       full_cmd,
                                       NUM_VARS,
                                       NUM_INDS);

    uint32_t i;
    for (i=0; i<NUM_VARS; ++i) 
        write_vid(v0, (i+1));

    destroy_vid_file(v0);

    struct vid_file *v1 = open_vid_file("test_vid");
    load_vid_data(v1);

    for (i=0; i<NUM_VARS; ++i)
        TEST_ASSERT_EQUAL(i+1, v1->vids[i]);

    destroy_vid_file(v1);

    remove("test_vid");
}
//}}}

void test_off_file(void)
{
}
