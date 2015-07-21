#include <limits.h>
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

//void test_check_field_name(void)
void test_check_field_name(void)
{
    // cannot start with a number
    TEST_ASSERT_EQUAL(0, check_field_name("123A"));

    // Test the invlid ranges
    int i;
    char test_str[3] = "A B";
    for (i = 32; i < 48; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 58; i < 65; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 91; i < 96; ++i) {
        if (i == '_')
            continue;
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 173; i < 128; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    TEST_ASSERT_EQUAL(-1, check_field_name("A"));
    TEST_ASSERT_EQUAL(-1, check_field_name("A234"));
    TEST_ASSERT_EQUAL(-1, check_field_name("A_234"));
}
//}}}

//{{{ void test_is_int(void)
void test_is_int(void)
{
    int v;
    TEST_ASSERT_EQUAL(0, is_int("134a", &v));
    TEST_ASSERT_EQUAL(0, is_int("134a", &v));

    TEST_ASSERT_EQUAL(0, is_int("134.3", &v));

    TEST_ASSERT_EQUAL(0, is_int(" 13 ", &v));

    // overflow
    TEST_ASSERT_EQUAL(0, is_int("928374238742387429340234", &v));
    TEST_ASSERT_EQUAL(0, is_int("2147483648", &v));


    TEST_ASSERT_EQUAL(1, is_int("13", &v));
    TEST_ASSERT_EQUAL(13, v);

    TEST_ASSERT_EQUAL(1, is_int(" 14", &v));
    TEST_ASSERT_EQUAL(14, v);

    TEST_ASSERT_EQUAL(1, is_int("2147483647", &v));
    TEST_ASSERT_EQUAL(INT_MAX, v);
}
//}}}
