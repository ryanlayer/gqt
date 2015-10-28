#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "unity.h"
#include "wahbm.h"
#include "parse_q.h"
#include "genotq.h"
#include "vid.h"
#include "bim.h"
#include "ubin.h"
#include "bcf.h"

char *BCF_FILE = "../data/diff_gts.bcf";
uint32_t NUM_INDS = 10;
uint32_t NUM_VARS = 43;

void setUp(void) { }
void tearDown(void) { }

//{{{void test_init_bcf_file_remote(void)
void test_init_bcf_file_remote(void)
{
    char *remote_file = "http://s3-us-west-2.amazonaws.com/gqt-data/test/10.1e4.var.bcf";
    struct bcf_file rbf = init_bcf_file(remote_file);
}
//}}}

//{{{void test_open_bim_file(void)
void test_open_bim_file(void)
{
    struct bim_file *b1 = open_bim_file("http://s3-us-west-2.amazonaws.com/gqt-data/test/10.1e4.var.bcf.bim");

    TEST_ASSERT_EQUAL('G', b1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', b1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', b1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('b', b1->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), b1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), b1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), b1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), b1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, b1->gqt_header->magic);

    TEST_ASSERT_EQUAL(43, b1->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(10, b1->gqt_header->num_samples);
#if 0
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
#endif

}
//}}}

//{{{void test_open_off_file(void)
void test_open_off_file(void)
{
    struct off_file *o1 = open_off_file("http://s3-us-west-2.amazonaws.com/gqt-data/test/10.1e4.var.bcf.off");

    TEST_ASSERT_EQUAL('G', o1->gqt_header->marker[0]);
    TEST_ASSERT_EQUAL('Q', o1->gqt_header->marker[1]);
    TEST_ASSERT_EQUAL('T', o1->gqt_header->marker[2]);
    TEST_ASSERT_EQUAL('o', o1->gqt_header->type);
    TEST_ASSERT_EQUAL(atoi(MAJOR_VERSION), o1->gqt_header->major);
    TEST_ASSERT_EQUAL(atoi(MINOR_VERSION), o1->gqt_header->minor);
    TEST_ASSERT_EQUAL(atoi(REVISION_VERSION), o1->gqt_header->revision);
    TEST_ASSERT_EQUAL(atoi(BUILD_VERSION), o1->gqt_header->build);
    TEST_ASSERT_EQUAL(0x11223344, o1->gqt_header->magic);

    TEST_ASSERT_EQUAL(43, o1->gqt_header->num_variants);
    TEST_ASSERT_EQUAL(10, o1->gqt_header->num_samples);
#if 0
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
#endif

    //remove(file_name);
}
//}}}
