#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "unity.h"
#include "parse_q.h"
#include "genotq.h"
#include "vid.h"

char *BCF_FILE = "../data/diff_gts.bcf";
uint32_t NUM_INDS = 10;
uint32_t NUM_VARS = 43;

void setUp(void) { }
void tearDown(void) { }

//void test_check_field_name(void)
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

void test_write_lod_vid_file(void)
{
    char *full_cmd = "gqt convert bcf -i bcf";
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
