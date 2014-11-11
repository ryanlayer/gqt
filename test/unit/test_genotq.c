#include <stdio.h>
#include <stdlib.h>
#include "parse_q.h"
#include "genotq.h"
#include "unity.h"
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include "timer.h"
#include "pthread_pool.h"
#include <string.h>

//{{{ SETUP
//{{{ char* itoa(int value, char* str, int radix) {
/* The itoa code is in the public domain */
char* itoa(int value, char* str, int radix) {
    static char dig[] =
        "0123456789"
        "abcdefghijklmnopqrstuvwxyz";
    int n = 0, neg = 0;
    unsigned int v;
    char* p, *q;
    char c;
    if (radix == 10 && value < 0) {
        value = -value;
        neg = 1;
    }
    v = value;
    do {
        str[n++] = dig[v%radix];
        v /= radix;
    } while (v);
    if (neg)
        str[n++] = '-';
    str[n] = '\0';
    for (p = str, q = p + (n-1); p < q; ++p, --q)
        c = *p, *p = *q, *q = c;
    return str;
}
//}}}

//{{{ int clear_list (struct wah_ll *A_head)
int clear_list (struct wah_ll *A_head)
{
    int c = 0;
    struct wah_ll *A_curr = A_head;
    while (A_curr != NULL) {
        struct wah_ll *A_tmp = A_curr->next;
        free(A_curr);
        A_curr = A_tmp;
        c += 1;
    }

    return c;
}
//}}}

//{{{ int clear_16list (struct wah_ll *A_head)
int clear_16list (struct wah16_ll *A_head)
{
    int c = 0;
    struct wah16_ll *A_curr = A_head;
    while (A_curr != NULL) {
        struct wah16_ll *A_tmp = A_curr->next;
        free(A_curr);
        A_curr = A_tmp;
        c += 1;
    }

    return c;
}
//}}}


void setUp(void)
{
}

void tearDown(void)
{
}
//}}}

//{{{void test_bin_char_to_int(void)
void test_bin_char_to_int(void)
{
    /*
     * 01000000000000000000000000000000 -> 0x40000000
     * 10000000000000000000000000000011 -> 0x80000003
     * 00011000000000000000000000000000 -> 0x18000000
     * 00000100000000000000000000000000 -> 0x04000000
     */

    char *c1 = "01000000000000000000000000000000";
    TEST_ASSERT_EQUAL(0x40000000, bin_char_to_int(c1));
    char *c2 = "10000000000000000000000000000011";
    TEST_ASSERT_EQUAL(0x80000003, bin_char_to_int(c2));
    char *c3 = "00011000000000000000000000000000";
    TEST_ASSERT_EQUAL(0x18000000, bin_char_to_int(c3));
    char *c4 = "00000100000000000000000000000000";
    TEST_ASSERT_EQUAL(0x04000000, bin_char_to_int(c4));

    TEST_ASSERT_EQUAL(0x40000000, 
                      bin_char_to_int("01000000000000000000000000000000"));
}
//}}}

//{{{void test_init_plt_file_num_fields_and_num_records(void)
void test_init_plt_file_num_fields_and_num_records(void)
{
    struct plt_file pltf_ind = init_plt_file("../data/10.1e4.ind.txt");
    TEST_ASSERT_EQUAL(43, pltf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, pltf_ind.num_records);
    fclose(pltf_ind.file);

    struct plt_file pltf_var = init_plt_file("../data/10.1e4.var.txt");
    TEST_ASSERT_EQUAL(10, pltf_var.num_fields);
    TEST_ASSERT_EQUAL(43, pltf_var.num_records);
    fclose(pltf_var.file);
}
//}}}

//{{{void test_pack_int(void)
void test_pack_int(void)
{
    int q1[16] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    unsigned int r =  pack_2_bit_ints(q1, 16);
    
    int q2[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r =  pack_2_bit_ints(q2, 16);
    TEST_ASSERT_EQUAL(1, r);

    int q3[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1};
    r =  pack_2_bit_ints(q3, 16);
    TEST_ASSERT_EQUAL(5, r);

    int q4[16] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r =  pack_2_bit_ints(q4, 16);
    TEST_ASSERT_EQUAL(pow(2,30)+1, r);

    int q5[16] = {2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3};
    r =  pack_2_bit_ints(q5, 16);
    TEST_ASSERT_EQUAL(pow(2,31)+3, r);

    int q6[16] = {2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r =  pack_2_bit_ints(q6, 16);
    TEST_ASSERT_EQUAL(pow(2,31)+pow(2,29)+1, r);
}
//}}}

//{{{ void test_unpack_int(void)
void test_unpack_int(void)
{
    int i, *r;

    int q2[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r = unpack_2_bit_ints(1);
    for (i = 0; i < 16; ++i) 
        TEST_ASSERT_EQUAL(q2[i], r[i]);
    free(r);

    int q3[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1};
    r = unpack_2_bit_ints(5);
    for (i = 0; i < 16; ++i)
        TEST_ASSERT_EQUAL(q3[i], r[i]);
    free(r);

    int q4[16] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r = unpack_2_bit_ints(pow(2,30)+1);
    for (i = 0; i < 16; ++i) 
        TEST_ASSERT_EQUAL(q4[i], r[i]);
    free(r);

    int q5[16] = {2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3};
    r = unpack_2_bit_ints(pow(2,31)+3);
    for (i = 0; i < 16; ++i) 
        TEST_ASSERT_EQUAL(q5[i], r[i]);
    free(r);

    int q6[16] = {2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    r = unpack_2_bit_ints(pow(2,31)+pow(2,29)+1);
    for (i = 0; i < 16; ++i) 
        TEST_ASSERT_EQUAL(q6[i], r[i]);
    free(r);
}
//}}}

//{{{ void test_pack_then_unpack_int(void)
void test_pack_then_unpack_int(void)
{
    int q[16], i, j, *r;

    for (j = 0; j < 100; ++j) {
        for (i = 0; i < 16; ++i)
            q[i] = rand()%4;
        
        unsigned int p = pack_2_bit_ints(q, 16);
        r = unpack_2_bit_ints(p);
        for (i = 0; i < 16; ++i) {
            TEST_ASSERT_EQUAL(q[i], r[i]);
        }
        free(r);
    }
}
//}}}

//{{{ void test_plt_line_to_packed_ints(void)
void test_plt_line_to_packed_ints(void)
{
    char *l1="2 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 2 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 0";
    //2 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 2 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 0
    //10 00 01 01 00 01 01 00 00 00 00 00 00 01 00 00 10 01 00 00 01 00 01 00 01 01 00 01 01 01 01 01 01 00 01 01 01 01 01 01 00 01 00
    //10 00 01 01 00 01 01 00 00 00 00 00 00 01 00 00
    //10 01 00 00 01 00 01 00 01 01 00 01 01 01 01 01
    //01 00 01 01 01 01 01 01 00 01 00
    //
    //10000101000101000000000000010000
    //10010000010001000101000101010101
    //0100010101010101000100
    //
    //10000101000101000000000000010000 -> 2232680464
    //10010000010001000101000101010101 -> 2420396373
    //01000101010101010001000000000000 -> 1163202560
   
    unsigned int *r;
    unsigned int len;

    len = plt_line_to_packed_ints(l1, 43, &r);
    TEST_ASSERT_EQUAL(len, 3);
    TEST_ASSERT_EQUAL(r[0], 2232680464);
    TEST_ASSERT_EQUAL(r[1], 2420396373);
    TEST_ASSERT_EQUAL(r[2], 1163202560);

    free(r);
    
    
    char *l2="1 0 0 0 0 0 1 1 1 0 1 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 1 0 1 1 1 1 1 1 0 1 0 1 1 0 1 0 0";
    //01000000000001010100010101000100 -> 1074087236
    //01000001000100010101010001010101 -> 1091654741
    //01010001000101000100000000000000 -> 1360281600
    
    len = plt_line_to_packed_ints(l2, 43, &r);
    TEST_ASSERT_EQUAL(len, 3);
    TEST_ASSERT_EQUAL(r[0], 1074087236);
    TEST_ASSERT_EQUAL(r[1], 1091654741);
    TEST_ASSERT_EQUAL(r[2], 1360281600);

    free(r);
}
//}}}

//{{{ void test_plot_to_ubin(void) 
void test_convert_file_plt_to_ubin(void) 
{

    //char *in_file_name="data/10.1e4.ind.txt";
    char *out_file_name="../data/.tmp";

    struct plt_file pltf_ind = init_plt_file("../data/10.1e4.ind.txt");
    convert_file_plt_to_ubin(pltf_ind, out_file_name);
    fclose(pltf_ind.file);
}
//}}}

//{{{ void test_init_ubin_file(void)
void test_init_ubin_file(void)
{
    struct plt_file pltf_ind = init_plt_file("../data/10.1e4.ind.txt");
    TEST_ASSERT_EQUAL(43, pltf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, pltf_ind.num_records);

    char *out_file_name_1="../data/.tmp1";
    convert_file_plt_to_ubin(pltf_ind, out_file_name_1);

    struct ubin_file uf_ind = init_ubin_file(out_file_name_1);
    TEST_ASSERT_EQUAL(43, uf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, uf_ind.num_records);

    fclose(uf_ind.file);


    struct plt_file pltf_var = init_plt_file("../data/10.1e4.var.txt");
    TEST_ASSERT_EQUAL(10, pltf_var.num_fields);
    TEST_ASSERT_EQUAL(43, pltf_var.num_records);

    char *out_file_name_2="../data/.tmp2";
    convert_file_plt_to_ubin(pltf_var, out_file_name_2);

    struct ubin_file uf_var = init_ubin_file(out_file_name_2);
    TEST_ASSERT_EQUAL(10, uf_var.num_fields);
    TEST_ASSERT_EQUAL(43, uf_var.num_records);

    fclose(uf_var.file);
    fclose(pltf_ind.file);
    fclose(pltf_var.file);
}
//}}} 

//{{{ void test_get_ubin_record(void)
void test_get_ubin_record(void)
{
    char *ubin_file_name="../data/10.1e4.ind.ubin";

    unsigned int A_0[43] = {2,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,2,1,0,0,
                            1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,
                            0,1,0};
    unsigned int A_1[43] = {1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,
                            0,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,
                            1,0,0};
    unsigned int A_2[43] = {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,0,2,
                            0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
                            0,1,0};
    unsigned int A_3[43] = {0,0,0,0,0,0,0,2,2,1,2,2,2,0,2,1,0,0,0,2,
                            0,2,0,1,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
                            0,0,0};
    unsigned int A_4[43] = {1,0,1,1,0,1,0,1,1,0,2,2,2,0,2,1,0,0,0,2,
                            0,2,0,2,0,0,1,0,0,2,0,0,2,2,0,0,0,0,0,0,
                            1,0,0};
    unsigned int A_5[43] = {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,1,2,
                            0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
                            0,0,0};
    unsigned int A_6[43] = {1,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,
                            0,0,0,0,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
                            0,0,0};
    unsigned int A_7[43] = {0,0,0,0,0,0,0,2,2,2,1,1,1,0,1,1,1,0,0,1,
                            0,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,
                            0,1,1};
    struct ubin_file uf = init_ubin_file(ubin_file_name);

    TEST_ASSERT_EQUAL(10, uf.num_records);
    TEST_ASSERT_EQUAL(43, uf.num_fields);

    unsigned int *ints, num_ints;
    num_ints = get_ubin_record(uf, 0, &ints);

    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);

    TEST_ASSERT_EQUAL(num_ints_per_record, num_ints);

    unsigned int i,j,two_bit,int_i = 0;
    for (i = 0; i < num_ints; ++i) {
        for (j = 0; j < 16; ++j) {
            two_bit = (ints[i] >> (30 - 2*j)) & 3;
            TEST_ASSERT_EQUAL(A_0[int_i], two_bit);
            int_i += 1;
            if (int_i == uf.num_fields)
                break;
        }
    }
    free(ints);

    num_ints = get_ubin_record(uf, 6, &ints);
    int_i = 0;
    for (i = 0; i < num_ints; ++i) {
        for (j = 0; j < 16; ++j) {
            two_bit = (ints[i] >> (30 - 2*j)) & 3;
            TEST_ASSERT_EQUAL(A_6[int_i], two_bit);
            int_i += 1;
            if (int_i == uf.num_fields)
                break;
        }
    }
    free(ints);

    num_ints = get_ubin_record(uf, 3, &ints);
    int_i = 0;
    for (i = 0; i < num_ints; ++i) {
        for (j = 0; j < 16; ++j) {
            two_bit = (ints[i] >> (30 - 2*j)) & 3;
            TEST_ASSERT_EQUAL(A_3[int_i], two_bit);
            int_i += 1;
            if (int_i == uf.num_fields)
                break;
        }
    }
    free(ints);


    fclose(uf.file);
}
//}}}

//{{{ void test_ints_to_rle(void)
void test_ints_to_rle(void)
{

    // 2^31 = 2147483648
    unsigned int I[5] = {2147483648,0,0,3,1};
    //10000000000000000000000000000000
    //00000000000000000000000000000000
    //00000000000000000000000000000000
    //00000000000000000000000000000011
    //00000000000000000000000000000001
    //1:1 0:125 1:2 0:31 1:1
    //1:1   -> 10000000000000000000000000000000 -> 2147483648
    //0:125 -> 00000000000000000000000001111100 -> 124
    //1:2   -> 10000000000000000000000000000001 -> 2147483649
    //0:31  -> 00000000000000000000000000011110 -> 30
    //1:1   -> 10000000000000000000000000000000 -> 2147483648
    unsigned int A[5] = {2147483648,124,2147483649,30,2147483648};

    unsigned int *O;
    int ints_to_rle_ints = ints_to_rle(I,5,&O);

    TEST_ASSERT_EQUAL(5, ints_to_rle_ints);

    int i;
    for (i= 0; i < 5; ++i)
        TEST_ASSERT_EQUAL(O[i], A[i]);

    free(O);

    unsigned int I2[6] = {0,0,0,0,0,0};

    ints_to_rle_ints = ints_to_rle(I2,6,&O);

    TEST_ASSERT_EQUAL(1, ints_to_rle_ints);
    TEST_ASSERT_EQUAL(O[0], 32*6-1);

    free(O);
}
//}}}

//{{{ void test_append_bit_to_active_word(void)
void test_append_bit_to_active_word(void)
{
    struct wah_active_word a;
    // 1101101 -> 7 bits, int = 109
    a.value = 109;
    a.nbits = 7; 

    int r = append_bit_to_active_word(&a, 1);

    TEST_ASSERT_EQUAL(0, r);
    // 11011011 -> 8 bits, int = 219
    TEST_ASSERT_EQUAL(8, a.nbits);
    TEST_ASSERT_EQUAL(219, a.value);
 
    r = append_bit_to_active_word(&a, 0);

    TEST_ASSERT_EQUAL(0, r);
    // 110110110 -> 9 bits, int = 438
    TEST_ASSERT_EQUAL(9, a.nbits);
    TEST_ASSERT_EQUAL(438, a.value);

    a.nbits = 31;
    a.value = 0;
    // Should covert to fill with fill bit set, the 
    // value bit set to zero and the #-> 32
    // 10....100000
    // 10000000000000000000000000100000 -> 2147483680
    r = append_bit_to_active_word(&a, 0);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(32, a.nbits);
    TEST_ASSERT_EQUAL(2147483680, a.value);

    r = append_bit_to_active_word(&a, 0);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(33, a.nbits);
    TEST_ASSERT_EQUAL(2147483681, a.value);

    a.nbits = 31;
    a.value = 2147483647;
    // Should covert to fill with fill bit set, the 
    // value bit set to one and the #-> 32
    // 11....100000
    // 11000000000000000000000000100000 -> 3221225504
    r = append_bit_to_active_word(&a, 1);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(32, a.nbits);
    TEST_ASSERT_EQUAL(3221225504, a.value);

    r = append_bit_to_active_word(&a, 1);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(33, a.nbits);
    TEST_ASSERT_EQUAL(3221225505, a.value);


    // try to add bits when the bit cannot be added
    // add a zero to a fill of ones
    a.nbits = 32;
    a.value = 3221225504;

    r = append_bit_to_active_word(&a, 0);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(32, a.nbits);
    TEST_ASSERT_EQUAL(3221225504, a.value);

    // add a one to a filee of zeros
    a.nbits = 33;
    a.value = 2147483681;

    r = append_bit_to_active_word(&a, 1);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(33, a.nbits);
    TEST_ASSERT_EQUAL(2147483681, a.value);
}
//}}}

//{{{void test_append_active_word(void)
void test_append_active_word(void)
{
    struct wah_ll *A_head = NULL,
                              *A_tail = NULL,
                              *A_curr;

    struct wah_active_word a;

    // Add to empty list
    a.nbits = 33;
    a.value = 2147483681;

    int r = append_active_word(&A_head,&A_tail,a);

    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(33, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(2147483681, A_tail->value.value);

    // Add a litteral with all zeros to to list containg a litteratal with all
    // a mix of zeros/ones
    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 31;
    //0101010101010101010101010101010 -> 715827882
    a.value = 715827882;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(715827882, A_tail->value.value);

    a.nbits = 31;
    a.value = 0;
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    //10000000000000000000000000111110 -> 
    TEST_ASSERT_EQUAL(0, A_tail->value.value);

    TEST_ASSERT_EQUAL(31, A_head->value.nbits);
    //10000000000000000000000000111110 -> 
    TEST_ASSERT_EQUAL(715827882, A_head->value.value);


    // Add a litteral with all zeros to to list containg a litteratal with all
    // zeros 
    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 31;
    a.value = 0;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(0, A_tail->value.value);

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(62, A_tail->value.nbits);
    //10000000000000000000000000000010 -> 0x80000002
    //TEST_ASSERT_EQUAL(2147483710, A_tail->value.value);
    TEST_ASSERT_EQUAL(0x80000002, A_tail->value.value);


    // Add a litteral with all zeros to to list containg a fill of 62  
    // zeros (2 words) 
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(93, A_tail->value.nbits);
    //10000000000000000000000000000011 -> 0x80000003
    TEST_ASSERT_EQUAL(0x80000003, A_tail->value.value);


    // Add a litteral with all ones
    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 31;
    // 1111111111111111111111111111111 -> 2147483647
    a.value = 2147483647;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(2147483647, A_tail->value.value);

    // Add a litteral with all ones where the list has a litteral with all ones
    a.nbits = 31;
    a.value = 2147483647;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(62, A_tail->value.nbits);
    // 11000000000000000000000000000010 -> 0xc0000002
    //TEST_ASSERT_EQUAL(3221225534, A_tail->value.value);
    TEST_ASSERT_EQUAL(0xc0000002, A_tail->value.value);

    // Add a litteral with all ones where the list has a fill of 31 ones
    // (1 word)
    a.nbits = 31;
    a.value = 2147483647;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(93, A_tail->value.nbits);
    // 11000000000000000000000000000011 -> 0xc0000003
    TEST_ASSERT_EQUAL(0xc0000003, A_tail->value.value);

    // Add a mixed litter to a list with a fill of ones
    a.nbits = 31;
    //0101010101010101010101010101010 -> 715827882
    a.value = 715827882;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(715827882, A_tail->value.value);


    // Add three mixed litterals
    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 31;
    a.value = 715827882;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(715827882, A_tail->value.value);

    a.nbits = 31;
    a.value = 715827882;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(715827882, A_tail->value.value);

    a.nbits = 31;
    a.value = 715827882;

    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(715827882, A_tail->value.value);

    A_curr = A_head;
    while (A_curr != NULL) {
        TEST_ASSERT_EQUAL(31, A_curr->value.nbits);
        TEST_ASSERT_EQUAL(715827882, A_curr->value.value);
        A_curr = A_curr->next;
    }
}
//}}}

//{{{ void test_map_from_32_bits_to_31_bits(void)
void test_map_from_32_bits_to_31_bits(void)
{

    // 2^31 = 2147483648
    unsigned int I[5] = {2147483648,0,0,3,1};
    /*
     * |-32---------------------------|
     * 10000000000000000000000000000000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000011
     * 00000000000000000000000000000001
     *
     * Regroup in to 31-bit groups
     *
     * |-31--------------------------|
     * 1000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0011000000000000000000000000000
     * 00001
     *
     * Pad the left size with zeros so each is 32-bits,
     * get the int values
     * |-32---------------------------|
     * 01000000000000000000000000000000 -> 1073741824
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00011000000000000000000000000000 -> 402653184
     * 00000100000000000000000000000000 -> 67108864
     */
    unsigned int A[6] = {1073741824, 0, 0, 0, 402653184, 67108864};
    unsigned int *O;
    unsigned int num_31_groups = map_from_32_bits_to_31_bits(I,5,160,&O);

    TEST_ASSERT_EQUAL(6, num_31_groups);

    int i;
    for (i= 0; i < 5; ++i)
        TEST_ASSERT_EQUAL(O[i], A[i]);

    free(O);

    /*
    * 
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
    
    num_31_groups = map_from_32_bits_to_31_bits(I2,6,192,&O);

    TEST_ASSERT_EQUAL(7, num_31_groups);

    for (i= 0; i < 7; ++i)
        TEST_ASSERT_EQUAL(O[i], A2[i]);

    free(O);
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

//{{{ void test_wah_run_decode(void)
void test_wah_run_decode(void)
{
    unsigned int I[5] = {2147483648,0,0,3,1};
    unsigned int A[4] = {0x40000000, 0x80000003, 0x18000000, 0x04000000};

    unsigned int *O;
    unsigned int wah_size = ints_to_wah(I,5,160,&O);

    TEST_ASSERT_EQUAL(4, wah_size);

    int i;
    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A[i]);

    struct wah_run r = init_wah_run(O, wah_size);

    TEST_ASSERT_EQUAL(0, r.word_i);
    TEST_ASSERT_EQUAL(wah_size, r.len);
    TEST_ASSERT_EQUAL(0, r.fill);
    TEST_ASSERT_EQUAL(0, r.num_words);
    TEST_ASSERT_EQUAL(0, r.is_fill);

    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(r.is_fill, 0);

    r.word_i += 1;
    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(1, r.is_fill);
    TEST_ASSERT_EQUAL(3, r.num_words);
    TEST_ASSERT_EQUAL(0, r.fill_bit);

    free(O);

    unsigned int I2[6] = {
            bin_char_to_int("10000000000000000000000000000000"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("00000000000000000000000000000011"),
            bin_char_to_int("00000000000000000000000000000001")
        };
 
    unsigned int A2[6] = {
            bin_char_to_int("01000000000000000000000000000000"),
            bin_char_to_int("00111111111111111111111111111111"),
            bin_char_to_int("11000000000000000000000000000010"),
            bin_char_to_int("01111000000000000000000000000000"),
            bin_char_to_int("00001100000000000000000000000000"),
            bin_char_to_int("00000010000000000000000000000000")
        };

    wah_size = ints_to_wah(I2,6,192,&O);

    TEST_ASSERT_EQUAL(6, wah_size);

    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A2[i]);

    r = init_wah_run(O, wah_size);

    TEST_ASSERT_EQUAL(0, r.word_i);
    TEST_ASSERT_EQUAL(wah_size, r.len);
    TEST_ASSERT_EQUAL(0, r.fill);
    TEST_ASSERT_EQUAL(0, r.num_words);
    TEST_ASSERT_EQUAL(0, r.is_fill);

    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(0, r.is_fill);

    r.word_i += 1;
    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(0, r.is_fill);
    TEST_ASSERT_EQUAL(1, r.num_words);

    r.word_i += 1;
    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(1, r.is_fill);
    TEST_ASSERT_EQUAL(2, r.num_words);
    TEST_ASSERT_EQUAL(1, r.fill_bit);

    r.word_i += 1;
    wah_run_decode(&r);
    TEST_ASSERT_EQUAL(0, r.is_fill);
    TEST_ASSERT_EQUAL(1, r.num_words);
}
//}}}

//{{{ void test_append_fill_word(void)
void test_append_fill_word(void)
{
    struct wah_ll *A_head = NULL,
                              *A_tail = NULL,
                              *A_curr;

    struct wah_active_word a;
    int r;


    //{{{ append 1-fill size 1 to empty list
    r = append_fill_word(&A_head,&A_tail,1,1);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1111111111111111111111111111111"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}

    //{{{ append 0-fill size 1 to empty list
    r = append_fill_word(&A_head,&A_tail,0,1);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("0000000000000000000000000000000"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}

    //{{{ append 1-fill size 2 to empty list
    r = append_fill_word(&A_head,&A_tail,1,2);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("11000000000000000000000000000010"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}

    //{{{ append 0-fill size 2 to empty list
    r = append_fill_word(&A_head,&A_tail,0,2);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("10000000000000000000000000000010"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}
   
    //{{{ append 1-fill size 1 to a list with mixed litterals
    a.nbits = 31;
    a.value = bin_char_to_int("1010100011101010101111101010111");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1010100011101010101111101010111"),
                      A_tail->value.value);


    r = append_fill_word(&A_head,&A_tail,1,1);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1111111111111111111111111111111"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(2, r);
    //}}}
    
    //{{{ append 1-fill size 2 to a list with mixed litterals
    a.nbits = 31;
    a.value = bin_char_to_int("1010100011101010101111101010111");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1010100011101010101111101010111"),
                      A_tail->value.value);


    r = append_fill_word(&A_head,&A_tail,1,2);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("11000000000000000000000000000010"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(2, r);
    //}}}

    //{{{ append 0-fill size 1 to a list with mixed litterals
    a.nbits = 31;
    a.value = bin_char_to_int("1010100011101010101111101010111");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1010100011101010101111101010111"),
                      A_tail->value.value);


    r = append_fill_word(&A_head,&A_tail,0,1);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("0000000000000000000000000000000"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(2, r);
    //}}}
    
    //{{{ append 0-fill size 2 to a list with mixed litterals
    a.nbits = 31;
    a.value = bin_char_to_int("1010100011101010101111101010111");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1010100011101010101111101010111"),
                      A_tail->value.value);


    r = append_fill_word(&A_head,&A_tail,0,2);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("10000000000000000000000000000010"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(2, r);
    //}}}

    //{{{ append 1-fill size 2 then 5  to a fill of 1s
    a.nbits = 31;
    a.value = bin_char_to_int("1111111111111111111111111111111");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("1111111111111111111111111111111"),
                      A_tail->value.value);
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("11000000000000000000000000000010"),
                      A_tail->value.value);

    r = append_fill_word(&A_head,&A_tail,1,2);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("11000000000000000000000000000100"),
                      A_tail->value.value);
 
    r = append_fill_word(&A_head,&A_tail,1,5);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("11000000000000000000000000001001"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}

    //{{{ append 0-fill size 2 then 5 to a fill of 0s
    a.nbits = 31;
    a.value = bin_char_to_int("0");
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(31, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(bin_char_to_int("0"),
                      A_tail->value.value);
    r = append_active_word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("10000000000000000000000000000010"),
                      A_tail->value.value);

    r = append_fill_word(&A_head,&A_tail,0,2);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("10000000000000000000000000000100"),
                      A_tail->value.value);
 
    r = append_fill_word(&A_head,&A_tail,0,5);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(bin_char_to_int("10000000000000000000000000001001"),
                      A_tail->value.value);

    r = clear_list(A_head);
    A_head = NULL;
    A_tail = NULL;
    TEST_ASSERT_EQUAL(1, r);
    //}}}
}
//}}}

//{{{ void test_wah_or(void)
void test_wah_or(void)
{
    /*
     * X
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1110100000000010101010000000000
     * 0000010000000000000000010101010
     * 00000
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * WAH
     * 536870912
     * 3221225474
     * 1946244096
     * 33554602
     * 0
     *
     * Y
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 01011
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     * WAH
     * 536870912
     * 3221225474
     * 2147483650
     * 738197504
     *
     *
     * Result:
     * 00100000000000000000000000000000
     * 11000000000000000000000000000010
     * 01110100000000010101010000000000
     * 00000010000000000000000010101010
     * 00101100000000000000000000000000
     */

    unsigned int A[5] = 
        { bin_char_to_int("00100000000000000000000000000000"),
          bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("01110100000000010101010000000000"),
          bin_char_to_int("00000010000000000000000010101010"),
          bin_char_to_int("00101100000000000000000000000000")
        };
 

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };

    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_or(&r_X, &r_Y, &Z);

    int i;
    for (i = 0; i < Z_len; ++i)
        TEST_ASSERT_EQUAL(A[i], Z[i]);

}
//}}}

//{{{ void test_wah_and(void)
void test_wah_and(void)
{
    /*
     * X
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * Y
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     *
     * Result:
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     */

    unsigned int A[5] = 
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")
        };
 

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };

    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_and(&r_X, &r_Y, &Z);

    unsigned int *ints;
    unsigned int ints_len = wah_to_ints(Z, Z_len, &ints);

    int i;
    for (i = 0; i < 5; ++i) 
        TEST_ASSERT_EQUAL(A[i], ints[i]);

}
//}}}

//{{{ void test_wah_to_ints(void)
void test_wah_to_ints(void)
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

    unsigned int *WAH;
    unsigned int wah_size = ints_to_wah(I,5,160,&WAH);

    unsigned int *INTS;
    unsigned int ints_size = wah_to_ints(WAH,wah_size,&INTS);

    unsigned int i;
    for (i = 0; i < 5; ++i)
        TEST_ASSERT_EQUAL(I[i], INTS[i]);

    free(INTS);
    free(WAH);

    unsigned int I2[6] = {1858530986,
                          2937432745,
                          1579791957,
                          3587518165,
                          1518381141,
                          1430257664};

    wah_size = ints_to_wah(I2,6,192,&WAH);
    ints_size = wah_to_ints(WAH,wah_size,&INTS);

    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(I2[i], INTS[i]);

    free(INTS);
    free(WAH);

    unsigned int I3[6] = {
            bin_char_to_int("10000000000000000000000000000000"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("11111111111111111111111111111111"),
            bin_char_to_int("00000000000000000000000000000011"),
            bin_char_to_int("00000000000000000000000000000001")
        };

    wah_size = ints_to_wah(I3,6,192,&WAH);
    ints_size = wah_to_ints(WAH,wah_size,&INTS);

    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(I3[i], INTS[i]);

    free(INTS);
    free(WAH);
}
//}}}

//{{{ void test_igned int *ints;
void test_ubin_to_bitmap(void)
{
    /*
     *  0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3
     *  0 0 1 1 1 1 2 2 3 3 0 0 1 1 1 2
     *  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     *  3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     *
     * 00011011000110110001101100011011
     * 00000101010110101111000001010110
     * 00000000000000000000000000000000
     * 11111111111111111111111111111111
     *
     * 0:
     * 10001000100010001100000000110000
     * 11111111111111110000000000000000
     * 1:
     * 01000100010001000011110000001110
     * 00000000000000000000000000000000
     * 2:
     * 00100010001000100000001100000001
     * 00000000000000000000000000000000
     * 3:
     * 00010001000100010000000011000000
     * 00000000000000001111111111111111
     *
     */

    unsigned int U1[4] = {
            bin_char_to_int("00011011000110110001101100011011"),
            bin_char_to_int("00000101010110101111000001010110"),
            bin_char_to_int("00000000000000000000000000000000"),
            bin_char_to_int("11111111111111111111111111111111")
    };

    unsigned int A1[8] = {
            bin_char_to_int("10001000100010001100000000110000"),
            bin_char_to_int("11111111111111110000000000000000"),
            bin_char_to_int("01000100010001000011110000001110"),
            bin_char_to_int("00000000000000000000000000000000"),
            bin_char_to_int("00100010001000100000001100000001"),
            bin_char_to_int("00000000000000000000000000000000"),
            bin_char_to_int("00010001000100010000000011000000"),
            bin_char_to_int("00000000000000001111111111111111")
    };

    unsigned int *B;

    unsigned int B_len = ubin_to_bitmap(U1,4,128,&B);

    TEST_ASSERT_EQUAL(8, B_len);
    unsigned int i;

    for (i = 0; i < 8; ++i) 
        TEST_ASSERT_EQUAL(A1[i], B[i]);

    free(B);

    unsigned int U2[3] = {
            bin_char_to_int("00011011000110110001101100011011"),
            bin_char_to_int("00000101010110101111000001010110"),
            bin_char_to_int("00000000000000000000000000000000")
    };

    unsigned int A2[8] = {
            // 0
            bin_char_to_int("10001000100010001100000000110000"),
            bin_char_to_int("11111111111111110000000000000000"),
            // 1
            bin_char_to_int("01000100010001000011110000001110"),
            bin_char_to_int("00000000000000000000000000000000"),
            // 2
            bin_char_to_int("00100010001000100000001100000001"),
            bin_char_to_int("00000000000000000000000000000000"),
            // 3
            bin_char_to_int("00010001000100010000000011000000"),
            bin_char_to_int("00000000000000000000000000000000")
    };

    B_len = ubin_to_bitmap(U2,3,96,&B);

    TEST_ASSERT_EQUAL(8, B_len);

    for (i = 0; i < 8; ++i) 
        TEST_ASSERT_EQUAL(A2[i], B[i]);


}
//}}}

//{{{void test_ubin_to_wah(void)
void test_ubin_to_wah(void)
{
    /*
     * int
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     *
     * ubin
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     *
     * wah
     * |---31------------------------|
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000101010101010101010101010101
     * 0101010101010101010101010101010
     * 1010101010101010101010101010101
     * 0101011010101010101010101010101
     * 0101010101010101010101010101010
     * 1010101010101010101010101010101
     * 0101010101111111111111111111111
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1111111111110000000000000000000
     *             ^-padding
     * 10000000000000000000000000000011
     * 00000101010101010101010101010101
     * 00101010101010101010101010101010
     * 01010101010101010101010101010101
     * 00101011010101010101010101010101
     * 00101010101010101010101010101010
     * 01010101010101010101010101010101
     * 00101010101111111111111111111111
     * 11000000000000000000000000000010
     * 01111111111110000000000000000000
     */
    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";
    unsigned int A[10] = {
            bin_char_to_int("10000000000000000000000000000011"),
            bin_char_to_int("00000101010101010101010101010101"),
            bin_char_to_int("00101010101010101010101010101010"),
            bin_char_to_int("01010101010101010101010101010101"),
            bin_char_to_int("00101011010101010101010101010101"),
            bin_char_to_int("00101010101010101010101010101010"),
            bin_char_to_int("01010101010101010101010101010101"),
            bin_char_to_int("00101010101111111111111111111111"),
            bin_char_to_int("11000000000000000000000000000010"),
            bin_char_to_int("01111111111110000000000000000000")
    };

    unsigned int *ubin;
    unsigned int ubin_len = plt_line_to_packed_ints(plt, 192, &ubin);

    TEST_ASSERT_EQUAL(12, ubin_len);

    unsigned int *wah;
    unsigned int wah_len = ints_to_wah(ubin, ubin_len, 192*2, &wah);

    TEST_ASSERT_EQUAL(10, wah_len);

    int i;
    for (i = 0; i < wah_len; ++i)
        TEST_ASSERT_EQUAL(A[i], wah[i]);


}
//}}}

//{{{ void test_ubin_to_bitmap_wah(void)
void test_ubin_to_bitmap_wah(void)
{
    /*
     * int
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     *
     * ubin
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     *
     * bit map
     * 0:
     * 11111111111111111111111111111111 -> 4294967295
     * 11111111111111110000000000000000 -> 4294901760
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 1: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000001111111111111111 -> 65535
     * 11111111111111111111111111111111 -> 4294967295
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 2: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 11111111111111111111111111111111 -> 4294967295
     * 11111111111111110000000000000000 -> 4294901760
     * 00000000000000000000000000000000 -> 0
     * 3: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000001111111111111111 -> 65535
     * 11111111111111111111111111111111 -> 4294967295
     *
     * wah
     * 0:
     * |--31-------------------------|
     * 1111111111111111111111111111111
     * 1111111111111111100000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     *       |-25 padding------------| 
     *
     * 01111111111111111111111111111111 -> 2147483647
     * 01111111111111111100000000000000 -> 2147467264
     * 10000000000000000000000000000101 -> 2147483653
     *
     * 1:
     * |--31-------------------------|
     * 0000000000000000000000000000000
     * 0000000000000000011111111111111
     * 1111111111111111111111111111111
     * 1110000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     *
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000011111111111111 -> 16383
     * 01111111111111111111111111111111 -> 2147483647
     * 01110000000000000000000000000000 -> 1879048192
     * 10000000000000000000000000000011 -> 2147483651
     * 2: 
     * |--31-------------------------|
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0001111111111111111111111111111
     * 1111111111111111111100000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     *
     * 10000000000000000000000000000011 -> 2147483651
     * 00001111111111111111111111111111 -> 268435455
     * 01111111111111111111100000000000 -> 2147481600
     * 10000000000000000000000000000010 -> 2147483650
     * 3: 
     * |--31-------------------------|
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 0000000000000000000011111111111
     * 1111111111111111111111111111111
     * 1111110000000000000000000000000
     *
     * 10000000000000000000000000000100 -> 2147483652
     * 00000000000000000000011111111111 -> 2047
     * 01111111111111111111111111111111 -> 2147483647
     * 01111110000000000000000000000000 -> 2113929216
     */

    unsigned int wah_A[16] = {
            //0
            2147483647,
            2147467264,
            2147483653,
            //1
            0,
            16383,
            2147483647,
            1879048192,
            2147483651,
            //2
            2147483651,
            268435455,
            2147481600,
            2147483650,
            //3
            2147483652,
            2047,
            2147483647,
            2113929216
    };

    unsigned int i;
    unsigned int wah_offsets_A[4] = {3,5,4,4};

    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";

    unsigned int *ubin;
    unsigned int ubin_len = plt_line_to_packed_ints(plt, 192, &ubin);

    TEST_ASSERT_EQUAL(12, ubin_len);

    unsigned int *wah;
    unsigned int *wah_offsets;
    unsigned int wah_len = ubin_to_bitmap_wah(ubin,
                                              ubin_len,
                                              192,
                                              &wah,
                                              &wah_offsets);

    TEST_ASSERT_EQUAL(16, wah_len);

    for (i = 0; i < 4; ++i)
        TEST_ASSERT_EQUAL(wah_offsets_A[i], wah_offsets[i]);

    for (i = 0; i < wah_len; ++i)
        TEST_ASSERT_EQUAL(wah_A[i], wah[i]);
}
//}}}

//{{{ void test_plt_to_bitmap_wah(void)
void test_plt_to_bitmap_wah(void)
{
    unsigned int wah_A[16] = {
            //0
            2147483647,
            2147467264,
            2147483653,
            //1
            0,
            16383,
            2147483647,
            1879048192,
            2147483651,
            //2
            2147483651,
            268435455,
            2147481600,
            2147483650,
            //3
            2147483652,
            2047,
            2147483647,
            2113929216
    };

    unsigned int i;
    unsigned int wah_offsets_A[4] = {3,5,4,4};

    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";

    unsigned int *wah;
    unsigned int *wah_offsets;
    unsigned int wah_len = plt_to_bitmap_wah(plt,
                                             192,
                                             &wah,
                                             &wah_offsets);

    TEST_ASSERT_EQUAL(16, wah_len);

    for (i = 0; i < 4; ++i)
        TEST_ASSERT_EQUAL(wah_offsets_A[i], wah_offsets[i]);

    for (i = 0; i < wah_len; ++i)
        TEST_ASSERT_EQUAL(wah_A[i], wah[i]);
}
//}}}

//{{{ void test_convert_file_by_name_ubin_to_wahbm(void)
void test_convert_file_by_name_ubin_to_wahbm(void)
{

    char *plt_file_name="../data/10.1e4.ind.txt";
    char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
    convert_file_by_name_ubin_to_wahbm(ubin_file_name, wah_file_name);

    struct wah_file wf = init_wahbm_file(wah_file_name);
    struct ubin_file uf = init_ubin_file(ubin_file_name);

    unsigned int test_record, test_bitmap;
    unsigned int *ints, num_ints;
    unsigned int *wah_bms[4], wah_sizes[4];
    unsigned int *wah_ints[4], wah_num_ints[4];
    unsigned int two_bit, bit, ubin_int_i, ubin_bit_i,
                 wah_int_i, wah_bit_i, field_i;
    for (test_record = 0; test_record < 8; ++test_record) {
        field_i = 0;
        num_ints = get_ubin_record(uf, test_record, &ints);

        for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap) {
            wah_sizes[test_bitmap] = get_wah_bitmap(wf,
                                                    test_record,
                                                    test_bitmap,
                                                    &(wah_bms[test_bitmap]));
            wah_num_ints[test_bitmap] = wah_to_ints(wah_bms[test_bitmap],
                                                    wah_sizes[test_bitmap],
                                                    &(wah_ints[test_bitmap]));

        }

        wah_int_i = 0;
        wah_bit_i = 0;

        for (ubin_int_i = 0; ubin_int_i < num_ints; ++ubin_int_i) {
            for (ubin_bit_i = 0; ubin_bit_i < 16; ++ubin_bit_i) {
                two_bit = (ints[ubin_int_i] >> (30-(ubin_bit_i*2))) & 3;

                for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap) {
                    bit = (wah_ints[test_bitmap][wah_int_i] >> 
                            (31-wah_bit_i)) & 1;

                    if (test_bitmap == two_bit)
                        TEST_ASSERT_EQUAL(1,bit);
                    else
                        TEST_ASSERT_EQUAL(0,bit);
                }

                wah_bit_i += 1;

                if (wah_bit_i == 32) {
                    wah_int_i += 1;
                    wah_bit_i = 0;
                }


                field_i += 1;
                if (field_i == wf.num_fields)
                    break;
            }
            if (field_i == wf.num_fields)
                break;
        }
        //fprintf(stderr, "\n");


        free(ints);
        for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap)
            free(wah_bms[test_bitmap]);
    }

    fclose(uf.file);
    fclose(wf.file);
    free(wf.record_offsets);
}
//}}}

//{{{ void test_convert_file_by_name_ubin_to_wah(void)
void test_convert_file_by_name_ubin_to_wah(void)
{

    char *plt_file_name="../data/10.1e4.ind.txt";
    char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wah";

    convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
    convert_file_by_name_ubin_to_wah(ubin_file_name, wah_file_name);

    struct wah_file wf = init_wah_file(wah_file_name);
    struct ubin_file uf = init_ubin_file(ubin_file_name);

    unsigned int field_count,
                 i,
                 j,
                 test_record,
                 num_ints,
                 num_wahs,
                 num_wah_ints;

    unsigned int *ints = NULL, *wahs = NULL, *wah_ints = NULL;

    for (test_record = 0; test_record < 8; ++test_record) {
        num_ints = get_ubin_record(uf, test_record, &ints);
        num_wahs = get_wah_record(wf, test_record, &wahs);
        num_wah_ints = wah_to_ints(wahs, num_wahs, &wah_ints);

        field_count = 0;
        for (i = 0; i < num_wah_ints; ++i) {
            for (j = 0; j < 16; ++j) {
                unsigned int wah_v = (wah_ints[i] >> (30 - 2*j)) & 3;
                unsigned int ubin_v = (ints[i] >> (30 - 2*j)) & 3;

                TEST_ASSERT_EQUAL(wah_v, ubin_v);
                field_count += 1;
                if (field_count == wf.num_fields )
                    break;
            }
            if (field_count == wf.num_fields ) 
                break;
        }

        free(ints);
        free(wahs);
        free(wah_ints);
        ints = NULL;
        wahs = NULL;
        wah_ints = NULL;
    }

    fclose(uf.file);
    fclose(wf.file);
    free(wf.record_offsets);
}
//}}}

//{{{ void test_convert_file_by_name_vcf_to_plt(void)
void test_convert_file_by_name_vcf_to_plt(void)
{

    char *vcf_file_name="../data/10.1e4.var.vcf";
    char *new_plt_file_name="../data/tmp.10.1e4.var.plt";
    char *orig_plt_file_name="../data/10.1e4.var.txt";

    convert_file_by_name_vcf_to_plt(vcf_file_name, 10, 43, new_plt_file_name);

    struct plt_file o_pf = init_plt_file(orig_plt_file_name);
    struct plt_file n_pf = init_plt_file(new_plt_file_name);

    TEST_ASSERT_EQUAL(o_pf.num_fields,n_pf.num_fields);
    TEST_ASSERT_EQUAL(o_pf.num_records,n_pf.num_records);

    unsigned int i, j, num_o_pf, num_n_pf;
    unsigned int *o_pfs = NULL, *n_pfs = NULL;

    for(i = 0; i < o_pf.num_records; ++i) {
        num_o_pf = get_plt_record(o_pf, i, &o_pfs);
        num_n_pf = get_plt_record(n_pf, i, &n_pfs);

        TEST_ASSERT_EQUAL(num_o_pf, num_n_pf);

        for(j = 0; j < num_o_pf; ++j) 
            TEST_ASSERT_EQUAL(o_pfs[j], n_pfs[j]);

        free(o_pfs);
        free(n_pfs);
        o_pfs = NULL;
        n_pfs = NULL;
    }
    

    fclose(o_pf.file);
    fclose(n_pf.file);

}
//}}}

//{{{ void test_init_wahbm_file(void)
void test_init_wahbm_file(void)
{
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct wah_file wf = init_wahbm_file(wah_file_name);

    TEST_ASSERT_EQUAL(10, wf.num_records);
    TEST_ASSERT_EQUAL(43, wf.num_fields);


    unsigned int A_record_offsets[40] = {
        2,  4,  6,  7,  9,  11, 12, 13, 15, 17, 
        19, 20, 22, 24, 26, 27, 29, 31, 33, 34,
        36, 38, 40, 41, 43, 45, 47, 48, 50, 52,
        54, 55, 57, 59, 61, 62, 64, 66, 68, 69
    };

    unsigned int i;
    for (i = 0; i < wf.num_records*4; ++i)
        TEST_ASSERT_EQUAL(A_record_offsets[i], wf.record_offsets[i]);

    fclose(wf.file);
    free(wf.record_offsets);
}
//}}}

//{{{void test_get_wah_bitmap(void)
void test_get_wah_bitmap(void)
{
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    unsigned int A[8][43] = {
        {2,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,2,1,0,0,
         1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,
         0,1,0},
        {1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,
         0,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,0,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,1,0},
        {0,0,0,0,0,0,0,2,2,1,2,2,2,0,2,1,0,0,0,2,
         0,2,0,1,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,1,1,0,1,0,1,1,0,2,2,2,0,2,1,0,0,0,2,
         0,2,0,2,0,0,1,0,0,2,0,0,2,2,0,0,0,0,0,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,1,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,
         0,0,0,0,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {0,0,0,0,0,0,0,2,2,2,1,1,1,0,1,1,1,0,0,1,
         0,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,
         0,1,1},
    };

    struct wah_file wf = init_wahbm_file(wah_file_name);

    TEST_ASSERT_EQUAL(10, wf.num_records);
    TEST_ASSERT_EQUAL(43, wf.num_fields);

    unsigned int *wah_bm, wah_size;
    unsigned int *ints, num_ints, bit;
    unsigned int test_record, test_bitmap, field_i, i, j;

    for (test_record = 0; test_record < 8; ++test_record) {
        for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap) {
            wah_size = get_wah_bitmap(wf,test_record,test_bitmap,&wah_bm);
            num_ints = wah_to_ints(wah_bm, wah_size, &ints);

            field_i = 0;
            for (i = 0; i < num_ints; ++i) {
                for (j = 32; j > 0; --j) {
                    bit = (ints[i] >> (j-1)) & 1;

                    if (A[test_record][field_i] == test_bitmap) 
                        TEST_ASSERT_EQUAL(1, bit);
                    else
                        TEST_ASSERT_EQUAL(0, bit);

                    field_i += 1;
                    if (field_i == wf.num_fields)
                        break;
                }
                if (field_i == wf.num_fields)
                    break;
            }

            free(wah_bm);
            free(ints);
        }
    }
    fclose(wf.file);
    free(wf.record_offsets);
}
//}}}

//{{{void test_get_wah_record(void)
void test_get_wah_record(void)
{

    unsigned int A[8][43] = {
        {2,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,2,1,0,0,
         1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,
         0,1,0},
        {1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,
         0,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,0,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,1,0},
        {0,0,0,0,0,0,0,2,2,1,2,2,2,0,2,1,0,0,0,2,
         0,2,0,1,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,1,1,0,1,0,1,1,0,2,2,2,0,2,1,0,0,0,2,
         0,2,0,2,0,0,1,0,0,2,0,0,2,2,0,0,0,0,0,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,1,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,
         0,0,0,0,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {0,0,0,0,0,0,0,2,2,2,1,1,1,0,1,1,1,0,0,1,
         0,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,
         0,1,1},
    };

    char *wah_file_name="../data/10.1e4.ind.wah";
    struct wah_file wf = init_wah_file(wah_file_name);

    TEST_ASSERT_EQUAL(10, wf.num_records);
    TEST_ASSERT_EQUAL(43, wf.num_fields);

    unsigned int field_count,
                 i,
                 j,
                 test_record,
                 num_wahs,
                 num_wah_ints;

    unsigned int *wahs = NULL, *wah_ints = NULL;

    for (test_record = 0; test_record < 8; ++test_record) {
        num_wahs = get_wah_record(wf, test_record, &wahs);
        num_wah_ints = wah_to_ints(wahs, num_wahs, &wah_ints);

        field_count = 0;
        for (i = 0; i < num_wah_ints; ++i) {
            for (j = 0; j < 16; ++j) {
                unsigned int wah_v = (wah_ints[i] >> (30 - 2*j)) & 3;

                TEST_ASSERT_EQUAL(A[test_record][field_count], wah_v);
                field_count += 1;
                if (field_count == wf.num_fields )
                    break;
            }
            if (field_count == wf.num_fields ) 
                break;
        }

        free(wahs);
        free(wah_ints);
        wahs = NULL;
        wah_ints = NULL;
    }

    fclose(wf.file);
    free(wf.record_offsets);
}
//}}}

//{{{void test_get_plt_record(void)
void test_get_plt_record(void)
{

    unsigned int A[8][43] = {
        {2,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,2,1,0,0,
         1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,
         0,1,0},
        {1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,
         0,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,0,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,1,0},
        {0,0,0,0,0,0,0,2,2,1,2,2,2,0,2,1,0,0,0,2,
         0,2,0,1,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,1,1,0,1,0,1,1,0,2,2,2,0,2,1,0,0,0,2,
         0,2,0,2,0,0,1,0,0,2,0,0,2,2,0,0,0,0,0,0,
         1,0,0},
        {0,0,0,0,0,0,0,2,2,0,2,2,2,0,2,0,0,0,1,2,
         0,2,0,2,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {1,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,
         0,0,0,0,0,0,0,0,0,2,0,0,2,2,0,0,0,0,0,0,
         0,0,0},
        {0,0,0,0,0,0,0,2,2,2,1,1,1,0,1,1,1,0,0,1,
         0,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,
         0,1,1},
    };

    char *plt_file_name="../data/10.1e4.ind.txt";
    struct plt_file pf = init_plt_file(plt_file_name);

    TEST_ASSERT_EQUAL(10, pf.num_records);
    TEST_ASSERT_EQUAL(43, pf.num_fields);

    unsigned int field_count,
                 i,
                 j,
                 test_record,
                 num_plts,
                 num_wah_ints;

    unsigned int *plts = NULL;

    for (test_record = 0; test_record < 8; ++test_record) {
        num_plts = get_plt_record(pf, test_record, &plts);

        field_count = 0;
        for (i = 0; i < num_plts; ++i) {
            for (j = 0; j < 16; ++j) {
                unsigned int plt_v = (plts[i] >> (30 - 2*j)) & 3;

                TEST_ASSERT_EQUAL(A[test_record][field_count], plt_v);
                field_count += 1;
                if (field_count == pf.num_fields )
                    break;
            }
            if (field_count == pf.num_fields ) 
                break;
        }

        free(plts);
        plts = NULL;
    }

    fclose(pf.file);
}
//}}}

//{{{ void test_gt_records_plt_ubin_wahbm(void)
void test_gt_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = gt_records_plt(pf, test_records, 4, 0, &pf_R);

    unsigned int *uf_R;
    unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 0, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_wahbm(wf, test_records, 4, 0, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

                
    /*
     * 1000001110111010100101011110111111010110100
     * 0000000220222020000202020000020022000000010
     * 0000000221222021000202010000020022000000000
     * 1011010110222021000202020010020022000000100
     *
     * 0000000110111010000101010000010011000000000
     *
     * 00000001101110100001010100000100 -> 28972292
     * 11000000000                      -> 1536
     */

    unsigned int A[2] = {28972292,1536};
    unsigned int shift[2] = {0,21};
    unsigned int i;
    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
}
//}}}

//{{{ void test_gte_records_plt_ubin_wahbm(void)
void test_gte_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    //char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    //struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = gte_records_plt(pf, test_records, 4, 1, &pf_R);

    //unsigned int *uf_R;
    //unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 0, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = gte_records_wahbm(wf, test_records, 4, 1, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);
                
    /*
     * 1000001110111010100101011110111111010110100
     * 0000000220222020000202020000020022000000010
     * 0000000221222021000202010000020022000000000
     * 1011010110222021000202020010020022000000100
     *
     * 0000000110111010000101010000010011000000000
     *
     * 00000001101110100001010100000100 -> 28972292
     * 11000000000                      -> 1536
     */

    unsigned int A[2] = {28972292,1536};
    unsigned int shift[2] = {0,21};
    unsigned int i;
    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
}
//}}}

//{{{ void test_lt_records_plt_ubin_wahbm(void)
void test_lt_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    //char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    //struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = lt_records_plt(pf, test_records, 4, 3, &pf_R);

    //unsigned int *uf_R;
    //unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 3, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = lt_records_wahbm(wf, test_records, 4, 3, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

    /*
    * genotype vectors
    * 1000001110111010100101011110111111010110100
    * 0000000220222020000202020000020022000000010
    * 0000000221222021000202010000020022000000000
    * 1011010110222021000202020010020022000000100               
    * 
    * lt(3)
    * 1111111111111111111111111111111111111111111
    *
    * 1111111111111111111111111111111111111111111000000000000000000000
    * |-32---------------------------||-32---------------------------|
    *                                            |-pad---------------|
    *
    * 11111111111111111111111111111111111111111110000000000000000000
    * |-31--------------------------||-31--------------------------|
    *                                            |-pad-------------|
    *
    * 1111111111111111111111111111111 -> 
    * 1111111111110000000000000000000
    * |-31--------------------------|
    *
    * WAH
    * 01111111111111111111111111111111 -> 2147483647
    * 01111111111110000000000000000000 -> 2146959360
    *
    * need to stitch the WAH together. second bit of second WAH block is tacked to first int
    *
    * 11111111111111111111111111111111 -> 4294967295
    * 11111111111000000000000000000000 -> 4292870144
    */

    unsigned int i;
    unsigned int A[2] = {4294967295,4292870144};
    unsigned int shift[2] = {0,0};
    for (i = 0; i < 2; ++i) {
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] & A[i]);
    }

    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i]);
}
//}}}

//{{{ void test_lte_records_plt_ubin_wahbm(void)
void test_lte_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    //char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    //struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = lte_records_plt(pf, test_records, 4, 2, &pf_R);

    //unsigned int *uf_R;
    //unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 3, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = lte_records_wahbm(wf, test_records, 4, 2, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

    /*
    * genotype vectors
    * 1000001110111010100101011110111111010110100
    * 0000000220222020000202020000020022000000010
    * 0000000221222021000202010000020022000000000
    * 1011010110222021000202020010020022000000100               
    * 
    * lte(2)
    * 1111111111111111111111111111111111111111111
    *
    * 1111111111111111111111111111111111111111111000000000000000000000
    * |-32---------------------------||-32---------------------------|
    *                                            |-pad---------------|
    *
    * 11111111111111111111111111111111111111111110000000000000000000
    * |-31--------------------------||-31--------------------------|
    *                                            |-pad-------------|
    *
    * 1111111111111111111111111111111 -> 
    * 1111111111110000000000000000000
    * |-31--------------------------|
    *
    * WAH
    * 01111111111111111111111111111111 -> 2147483647
    * 01111111111110000000000000000000 -> 2146959360
    *
    * need to stitch the WAH together. second bit of second WAH block is tacked to first int
    *
    * 11111111111111111111111111111111 -> 4294967295
    * 11111111111000000000000000000000 -> 4292870144
    */

    unsigned int i;
    unsigned int A[2] = {4294967295,4292870144};
    unsigned int shift[2] = {0,0};
    for (i = 0; i < 2; ++i) {
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] & A[i]);
    }

    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i]);
}
//}}}

//{{{ void test_eq_records_plt_ubin_wahbm(void)
void test_eq_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    //char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    //struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = eq_records_plt(pf, test_records, 4, 0, &pf_R);

    //unsigned int *uf_R;
    //unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 3, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = eq_records_wahbm(wf, test_records, 4, 0, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

    /*
    * genotype vectors
    * 1000001110111010100101011110111111010110100
    * 0000000220222020000202020000020022000000010
    * 0000000221222021000202010000020022000000000
    * 1011010110222021000202020010020022000000100               
    *
    * eq(0)
    * 0100100000000100011010100001000000101001001
    *
    * 01001000000001000110101000010000 -> 1208248848
    * 00101001001                      -> 329
    */

    unsigned int i;
    unsigned int A[2] = {1208248848,329};
    unsigned int shift[2] = {0,21};
    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
}
//}}}

//{{{ void test_ne_records_plt_ubin_wahbm(void)
void test_ne_records_plt_ubin_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    //char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    //struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = ne_records_plt(pf, test_records, 4, 1, &pf_R);

    //unsigned int *uf_R;
    //unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 3, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = ne_records_wahbm(wf, test_records, 4, 1, &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

    /*
    * genotype vectors
    * 1000001110111010100101011110111111010110100
    * 0000000220222020000202020000020022000000010
    * 0000000221222021000202010000020022000000000
    * 1011010110222021000202020010020022000000100               
    * 
    * ne(1)
    * 0100100000000100011010100001000000101001001
    *
    * 01001000000000000110101000010000 -> 1208248848
    * 00101001001                      -> 329
    */

    unsigned int i;
    unsigned int A[2] = {1208248848,329};
    unsigned int shift[2] = {0,21};
    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    //for (i = 0; i < 2; ++i)
    //    TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
}
//}}}

//{{{void test_padding_fix(void)
void test_padding_fix(void)
{
    /*
     * plt:
     * 2011011000000100210010101101111110111111010
     *
     * 2011011000000100
     * 2100101011011111
     * 10111111010
     *
     * int:
     * 10000101000101000000000000010000 -> 2232680464
     * 10010000010001000101000101010101 -> 2420396373
     * 01000101010101010001000000000000 -> 1163202560
     *                       |-pad----|
     *
     * bit maps:
     * 0100100111111011001101010010000001000000101
     * 0011011000000100010010101101111110111111010
     * 1000000000000000100000000000000000000000000
     * 0000000000000000000000000000000000000000000
     *
     * 01001001111110110011010100100000 01000000101000000000000000000000
     * 00110110000001000100101011011111 10111111010000000000000000000000
     * 10000000000000001000000000000000 00000000000000000000000000000000
     * 00000000000000000000000000000000 00000000000000000000000000000000
     * |-32---------------------------| |-32---------------------------|
     *                                             |--22 bit int pad---|
     *
     * 31-bit groups: 
     * 1:
     * 0100100111111011001101010010000
     * 0010000001010000000000000000000
     * 2:
     * 0011011000000100010010101101111
     * 1101111110100000000000000000000
     * 3:
     * 1000000000000000100000000000000
     * 0000000000000000000000000000000
     * 4:
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * |-31--------------------------|
     *
     * WAH
     * 1:
     * 00100100111111011001101010010000 -> 
     * 00010000001010000000000000000000
     * 2:
     * 00011011000000100010010101101111
     * 01101111110100000000000000000000
     * 3:
     * 01000000000000000100000000000000
     * 00000000000000000000000000000000
     * 4:
     * 10000000000000000000000000000010
     * 
     */

    unsigned int i,j;
    char *plt = "2 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 "
                "2 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 "
                "1 0 1 1 1 1 1 1 0 1 0";

    unsigned int int_bm_A[8] = {
            1241199904,
            1084227584,
            906250975,
            3208642560,
            2147516416,
            0,
            0,
            0};

    unsigned int wahs_A[7] = {
            bin_char_to_int("00100100111111011001101010010000"),
            bin_char_to_int("00010000001010000000000000000000"),
            bin_char_to_int("00011011000000100010010101101111"),
            bin_char_to_int("01101111110100000000000000000000"),
            bin_char_to_int("01000000000000000100000000000000"),
            bin_char_to_int("00000000000000000000000000000000"),
            bin_char_to_int("10000000000000000000000000000010")
    };
    
    unsigned int *ints;
    unsigned int int_len = plt_line_to_packed_ints(plt, 43, &ints);

    TEST_ASSERT_EQUAL(3, int_len);

    unsigned int *int_bm;
    unsigned int int_bm_len = ubin_to_bitmap(ints, int_len, 86, &int_bm);

    TEST_ASSERT_EQUAL(8, int_bm_len);

    for (i = 0; i < int_bm_len; ++i)
        TEST_ASSERT_EQUAL(int_bm_A[i], int_bm[i]);

    unsigned int bm_len = int_bm_len / 4;
    unsigned int wah_lens[4], *wahs[4];
    for (i = 0; i < 4; ++i)
        wah_lens[i] = ints_to_wah(&(int_bm[i*2]), bm_len, 43, &(wahs[i]));

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < wah_lens[i]; ++j) {
            TEST_ASSERT_EQUAL(wahs_A[i*2+j], wahs[i][j]);
        }
    }

}
//}}}

//{{{void test_get_wah_bitmap_in_place(void)
void test_get_wah_bitmap_in_place(void)
{
    unsigned int i,j,k;

    char *wah_file_name="../data/10.1e4.ind.wahbm";
    struct wah_file wf1 = init_wahbm_file(wah_file_name);
    struct wah_file wf2 = init_wahbm_file(wah_file_name);

    unsigned int *wah_bm1;

    unsigned int max_wah_size = (wf2.num_fields + 31 - 1)/ 31;
    unsigned int *wah_bm2 = (unsigned int *)
            malloc(sizeof(unsigned int)*max_wah_size);


    for (i = 0; i < 10; ++i) {
        for (j = 0; j < 4; ++j) {
            unsigned int wah_size1 = get_wah_bitmap(wf1,i,j,&wah_bm1);
            unsigned int wah_size2 = get_wah_bitmap_in_place(wf2,i,j,&wah_bm2);

            TEST_ASSERT_EQUAL(wah_size1, wah_size2);

            for (k = 0; k < wah_size2; ++k)
                TEST_ASSERT_EQUAL(wah_bm1[k],wah_bm2[k]);

            free(wah_bm1);
        }
    }

    free(wah_bm2);
}
//}}}

//{{{ void test_get_wah_bitmap_in_place_then_or(void)
void test_get_wah_bitmap_in_place_then_or(void)
{
    unsigned int i,j,k;

    char *wah_file_name="../data/10.1e4.ind.wahbm";
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int max_wah_size = (wf.num_fields + 31 - 1)/ 31;

    unsigned int *wah_ip_bm1 = (unsigned int *)
            malloc(sizeof(unsigned int)*max_wah_size);
    unsigned int *wah_ip_bm2 = (unsigned int *)
            malloc(sizeof(unsigned int)*max_wah_size);

    unsigned int wah_ip_size1 = get_wah_bitmap_in_place(wf,1,2,&wah_ip_bm1);
    unsigned int wah_ip_size2 = get_wah_bitmap_in_place(wf,2,2,&wah_ip_bm2);

    struct wah_run r_ip_bm1 = init_wah_run(wah_ip_bm1, wah_ip_size1);
    struct wah_run r_ip_bm2 = init_wah_run(wah_ip_bm2, wah_ip_size2);

    unsigned int *wah_bm1, *wah_bm2;

    unsigned int wah_size1 = get_wah_bitmap(wf,1,2,&wah_bm1);
    unsigned int wah_size2 = get_wah_bitmap(wf,2,2,&wah_bm2);

    struct wah_run r_bm1 = init_wah_run(wah_bm1, wah_size1);
    struct wah_run r_bm2 = init_wah_run(wah_bm2, wah_size2);


    unsigned int *r_ip;
    unsigned int r_ip_len = wah_or(&r_ip_bm1, &r_ip_bm2, &r_ip);
    
    unsigned int *r;
    unsigned int r_len = wah_or(&r_bm1, &r_bm2, &r);

    TEST_ASSERT_EQUAL(r_len, r_ip_len);

    for (i = 0; i < r_len; ++i) 
        TEST_ASSERT_EQUAL(r[i], r_ip[i]);
}
//}}}

//{{{ void test_map_from_32_bits_to_16_bits(void)
void test_map_from_32_bits_to_16_bits(void)
{
    /*
     * 00100100111111011001101010010000
     * 00010000001010000000000000000000
     * 00011011000000100010010101101111
     * 01101111110100000000000000000000
     * 01000000000000000100000000000000
     * 00000000000000000000000000000000
     * 10000000000000000000000000000010
     * |-32---------------------------|
     *
     * |-15----------|
     * 001001001111110
     * 110011010100100
     * 000001000000101
     * 000000000000000
     * 000000011011000
     * 000100010010101
     * 101111011011111
     * 101000000000000
     * 000000000100000
     * 000000000010000
     * 000000000000000
     * 000000000000000
     * 000000000000100
     * 000000000000000
     * 000000000000100
     *               ^pad     
     * 
     */
 

    unsigned int ints[7] = {
            bin_char_to_int("00100100111111011001101010010000"),
            bin_char_to_int("00010000001010000000000000000000"),
            bin_char_to_int("00011011000000100010010101101111"),
            bin_char_to_int("01101111110100000000000000000000"),
            bin_char_to_int("01000000000000000100000000000000"),
            bin_char_to_int("00000000000000000000000000000000"),
            bin_char_to_int("10000000000000000000000000000010")
    };

    unsigned int O_A[15] = {
        bin_char_to_int("001001001111110"),
        bin_char_to_int("110011010100100"),
        bin_char_to_int("000001000000101"),
        bin_char_to_int("000000000000000"),
        bin_char_to_int("000000011011000"),
        bin_char_to_int("000100010010101"),
        bin_char_to_int("101111011011111"),
        bin_char_to_int("101000000000000"),
        bin_char_to_int("000000000100000"),
        bin_char_to_int("000000000010000"),
        bin_char_to_int("000000000000000"),
        bin_char_to_int("000000000000000"),
        bin_char_to_int("000000000000100"),
        bin_char_to_int("000000000000000"),
        bin_char_to_int("000000000000100")
    };

    unsigned int i;
    uint16_t *O;
    unsigned int O_len = map_from_32_bits_to_15_bits(ints,
                                                     7,
                                                     7*32,
                                                     &O);
    TEST_ASSERT_EQUAL(15, O_len);

    for (i = 0; i < O_len; ++i) 
        TEST_ASSERT_EQUAL(O_A[i], O[i]);
}
//}}}

//{{{void test_append_active_16word(void)
void test_append_active_16word(void)
{
    struct wah16_ll *A_head = NULL,
                    *A_tail = NULL,
                    *A_curr;

    struct wah16_active_word a;

    // Add to empty list
    // 100101010101010 -> 19114
    a.nbits = 15;
    a.value = 19114;

    int r = append_active_16word(&A_head,&A_tail,a);

    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(19114, A_tail->value.value);

    // Add a litteral with all zeros to to list containg a litteratal with all
    // a mix of zeros/ones
    r = clear_16list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 15;
    //111000111000111 -> 29127
    a.value = 29127;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(29127, A_tail->value.value);

    a.nbits = 15;
    a.value = 0;
    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    //10000000000000000000000000111110 -> 
    TEST_ASSERT_EQUAL(0, A_tail->value.value);

    TEST_ASSERT_EQUAL(15, A_head->value.nbits);
    //10000000000000000000000000111110 -> 
    TEST_ASSERT_EQUAL(29127, A_head->value.value);


    // Add a litteral with all zeros to to list containg a litteratal with all
    // zeros 
    r = clear_16list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 15;
    a.value = 0;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(0, A_tail->value.value);

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //1000000000000010 -> 0x8002
    TEST_ASSERT_EQUAL(0x8002, A_tail->value.value);

    // Add a litteral with all zeros to to list containg a fill of 62  
    // zeros (2 words) 
    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(93, A_tail->value.nbits);
    //10000000000000000000000000000011 -> 0x80000003
    TEST_ASSERT_EQUAL(0x8003, A_tail->value.value);


    // Add a litteral with all ones
    r = clear_16list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 15;
    // 111111111111111 -> 0x7fff
    a.value = 0x7fff;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(0x7fff, A_tail->value.value);

    // Add a litteral with all ones where the list has a litteral with all ones
    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    //TEST_ASSERT_EQUAL(62, A_tail->value.nbits);
    // 1100000000000010 -> 0xc002
    TEST_ASSERT_EQUAL(0xc002, A_tail->value.value);

    // Add a litteral with all ones where the list has a fill of 31 ones
    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(0, r);
    TEST_ASSERT_EQUAL(0xc003, A_tail->value.value);

    // Add a mixed litter to a list with a fill of ones
    a.nbits = 15;
    //1010101010101010 -> 43690
    a.value = 43690;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(43690, A_tail->value.value);


    // Add three mixed litterals
    r = clear_16list(A_head);
    A_head = NULL;
    A_tail = NULL;

    a.nbits = 15;
    a.value = 4151;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(4151, A_tail->value.value);

    a.nbits = 15;
    a.value = 30001;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(30001, A_tail->value.value);

    a.nbits = 15;
    a.value = 15115;

    r = append_active_16word(&A_head,&A_tail,a);
    TEST_ASSERT_EQUAL(1, r);
    TEST_ASSERT_EQUAL(15, A_tail->value.nbits);
    TEST_ASSERT_EQUAL(15115, A_tail->value.value);
}
//}}}

//{{{ void test_ints_to_wah16(void)
void test_ints_to_wah16(void)
{
    /*
     * |-32---------------------------|
     * 10000000000000000000000000000000 -> 2147483648
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000011 -> 3
     * 00000000000000000000000000000001 -> 1
     *
     * Regroup into 15-bit groups
     *
     *
     * 100000000000000 
     * 000000000000000 +
     * 000000000000000 |
     * 000000000000000 |
     * 000000000000000 |
     * 000000000000000 |
     * 000000000000000 |
     * 000000000000000 +
     * 000000110000000
     * 000000000000000
     * 000000000100000
     * |-15----------|
     *           |pad|
     *
     * WAH:
     * 0100000000000000 -> 16384
     * 1000000000000111 -> 32775
     * 0000000110000000 -> 384
     * 0000000000000000 -> 0
     * 0000000000100000 -> 32
     */
    unsigned int I[5] = {2147483648,0,0,3,1};
    unsigned int A[5] = {16384,32775,384,0,32};

    uint16_t *O;
    unsigned int wah_size = ints_to_wah16(I,5,160,&O);

    TEST_ASSERT_EQUAL(5, wah_size);

    int i;
    for (i = 0; i < wah_size; ++i)
        TEST_ASSERT_EQUAL(O[i], A[i]);
}
//}}}

//{{{ void test_ubin_to_bitmap_wah16(void)
void test_ubin_to_bitmap_wah16(void)
{
    /*
     * int
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     * 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
     *
     * ubin
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 01010101010101010101010101010101
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 10101010101010101010101010101010
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     *
     * bit map
     * 0:
     * 11111111111111111111111111111111 -> 4294967295
     * 11111111111111110000000000000000 -> 4294901760
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 1: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000001111111111111111 -> 65535
     * 11111111111111111111111111111111 -> 4294967295
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 2: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 11111111111111111111111111111111 -> 4294967295
     * 11111111111111110000000000000000 -> 4294901760
     * 00000000000000000000000000000000 -> 0
     * 3: 
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000000000000000000000 -> 0
     * 00000000000000001111111111111111 -> 65535
     * 11111111111111111111111111111111 -> 4294967295
     *
     * 15 bit groups
     * 0:
     * |-15----------|
     * 111111111111111
     * 111111111111111
     * 111111111111111
     * 111000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     *             |p|
     * 1: 
     * |-15----------|
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000111111111111
     * 111111111111111
     * 111111111111111
     * 111111000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     *             |p|
     *
     * 2: 
     * |-15----------|
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000111111111
     * 111111111111111
     * 111111111111111
     * 111111111000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     *             |p|
     * 3: 
     * |-15----------|
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000000000
     * 000000000111111
     * 111111111111111
     * 111111111111111
     * 111111111111000
     *             |p|
     * 
     * WAH
     * 0:
     * 1100000000000011 -> 49155
     * 0111000000000000 -> 28672
     * 1000000000001001 -> 32777
     * 1:              
     * 1000000000000011 -> 32771
     * 0000111111111111 -> 4095
     * 1100000000000010 -> 49154
     * 0111111000000000 -> 32256
     * 1000000000000110 -> 32774
     * 2:              
     * 1000000000000110 -> 32774
     * 0000000111111111 -> 511
     * 1100000000000010 -> 49154
     * 0111111111000000 -> 32704
     * 1000000000000011 -> 32771
     * 3:              
     * 1000000000001001 -> 32777
     * 0000000000111111 -> 63
     * 1100000000000010 -> 49154
     * 0111111111111000 -> 32760
     */

    unsigned int wah_A[17] = {
            //0
            49155,
            28672,
            32777,
            //1
            32771,
            4095,
            49154,
            32256,
            32774,
            //2
            32774,
            511,
            49154,
            32704,
            32771,
            //3
            32777,
            63,
            49154,
            32760
    };

    unsigned int i;
    unsigned int wah_offsets_A[4] = {3,5,5,4};

    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";

    unsigned int *ubin;
    unsigned int ubin_len = plt_line_to_packed_ints(plt, 192, &ubin);

    TEST_ASSERT_EQUAL(12, ubin_len);

    uint16_t *wah;
    unsigned int *wah_offsets;
    unsigned int wah_len = ubin_to_bitmap_wah16(ubin,
                                              ubin_len,
                                              192,
                                              &wah,
                                              &wah_offsets);

    TEST_ASSERT_EQUAL(17, wah_len);

    for (i = 0; i < 4; ++i)
        TEST_ASSERT_EQUAL(wah_offsets_A[i], wah_offsets[i]);

    for (i = 0; i < wah_len; ++i)
        TEST_ASSERT_EQUAL(wah_A[i], wah[i]);

}
//}}}

//{{{ void test_wah_in_place_or(void)
void test_wah_in_place_or(void)
{
    /*
     * X
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1110100000000010101010000000000
     * 0000010000000000000000010101010
     * 00000
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * WAH
     * 536870912
     * 3221225474
     * 1946244096
     * 33554602
     * 0
     *
     * Y
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 01011
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     * WAH
     * 536870912
     * 3221225474
     * 2147483650
     * 738197504
     *
     *
     * Result:
     * 00100000000000000000000000000000
     * 11000000000000000000000000000010
     * 01110100000000010101010000000000
     * 00000010000000000000000010101010
     * 00101100000000000000000000000000
     */

    unsigned int A[5] = 
        { bin_char_to_int("00100000000000000000000000000000"),
          bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("01110100000000010101010000000000"),
          bin_char_to_int("00000010000000000000000010101010"),
          bin_char_to_int("00101100000000000000000000000000")
        };
 

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };

    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_or(&r_X, &r_Y, &Z);

    int i;
    for (i = 0; i < Z_len; ++i)
        TEST_ASSERT_EQUAL(A[i], Z[i]);

    
    unsigned int R[6] = {0,0,0,0,0,0};

    unsigned int r = wah_in_place_or(R, 6, w_Y, wah_size_Y);
    r = wah_in_place_or(R, 6, w_X, wah_size_X);

    unsigned int *ints_o, *ints_ip;
    unsigned int ints_ip_len = wah_to_ints(R,r,&ints_ip);
    unsigned int ints_o_len = wah_to_ints(Z,Z_len,&ints_o);

    TEST_ASSERT_EQUAL(ints_o_len, ints_ip_len);
    for (i = 0; i < ints_o_len; ++i)
        TEST_ASSERT_EQUAL(ints_o[i], ints_ip[i]);

}
//}}}

//{{{ void test_wah_compressed_in_place_or(void)
void test_wah_compressed_in_place_or(void)
{
    /*
     * X
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1110100000000010101010000000000
     * 0000010000000000000000010101010
     * 00000
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * WAH
     * 536870912
     * 3221225474
     * 1946244096
     * 33554602
     * 0
     *
     * Y
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 01011
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     * WAH
     * 536870912
     * 3221225474
     * 2147483650
     * 738197504
     *
     *
     * Result:
     * 00100000000000000000000000000000
     * 11000000000000000000000000000010
     * 01110100000000010101010000000000
     * 00000010000000000000000010101010
     * 00101100000000000000000000000000
     */

    unsigned int A[5] = 
        { bin_char_to_int("00100000000000000000000000000000"),
          bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("01110100000000010101010000000000"),
          bin_char_to_int("00000010000000000000000010101010"),
          bin_char_to_int("00101100000000000000000000000000")
        };
 

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };

    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_or(&r_X, &r_Y, &Z);

    unsigned int i;
    for (i = 0; i < Z_len; ++i)
        TEST_ASSERT_EQUAL(A[i], Z[i]);

    // init the first value to be a fill of zeros as long as the max size is
    unsigned int R[6] = {0x80000006,0,0,0,0,0};


    /*
     * - 00100000000000000000000000000000  - 10000000000000000000000000000110
     *   11000000000000000000000000000010    0
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    0
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     * - 11000000000000000000000000000010  - 10000000000000000000000000000101
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    0
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     * - 10000000000000000000000000000010    0
     *   00101100000000000000000000000000  - 10000000000000000000000000000011
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   10000000000000000000000000000010    0
     * - 00101100000000000000000000000000    10000000000000000000000000000010
     *                                       0
     *                                     - 00000000000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    10000000000000000000000000000010
     *                                       0
     *                                     - 00101100000000000000000000000000
     *
     *   00100000000000000000000000000000 -> 536870912
     *   11000000000000000000000000000010 -> 3221225474
     *   0                                -> 0
     *   10000000000000000000000000000010 -> 2147483650
     *   0                                -> 0
     *   00101100000000000000000000000000 -> 738197504
     */
    
    unsigned int A0[6] = {536870912,
                          3221225474,
                          0,
                          2147483650,
                          0,
                          738197504};

    unsigned int r = wah_compressed_in_place_or(R, 6, w_Y, wah_size_Y);

    TEST_ASSERT_EQUAL(6, r);
        
    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(A0[i], R[i]);

    r = wah_compressed_in_place_or(R, 6, w_X, wah_size_X);

    /*
     * - 00100000000000000000000000000000  - 00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     * 
     *   00100000000000000000000000000000    00100000000000000000000000000000
     * - 11000000000000000000000000000010  - 11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     * - 01110100000000010101010000000000    0
     *   00000010000000000000000010101010  - 10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     * - 00000010000000000000000010101010    01110100000000010101010000000000
     *   00000000000000000000000000000000  - 00000000000000000000000000000000
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    01110100000000010101010000000000
     * - 00000000000000000000000000000000    00000010000000000000000010101010
     *                                     - 00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    01110100000000010101010000000000
     *   00000000000000000000000000000000    00000010000000000000000010101010
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000 -> 536870912
     *   11000000000000000000000000000010 -> 3221225474
     *   0                                -> 0
     *   01110100000000010101010000000000 -> 1946244096
     *   00000010000000000000000010101010 -> 33554602
     *   00101100000000000000000000000000 -> 738197504
     */

    unsigned int A1[6] = {536870912,
                          3221225474,
                          0,
                          1946244096,
                          33554602,
                          738197504};

    TEST_ASSERT_EQUAL(6, r);
        
    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(A1[i], R[i]);


    unsigned int *ints_o, *ints_ip;
    unsigned int ints_ip_len = compressed_in_place_wah_to_ints(R,r,&ints_ip);
    unsigned int ints_o_len = wah_to_ints(Z,Z_len,&ints_o);

    TEST_ASSERT_EQUAL(ints_o_len, ints_ip_len);

    for (i = 0; i < ints_o_len; ++i)
        TEST_ASSERT_EQUAL(ints_o[i], ints_ip[i]);
}
//}}}

//{{{ void test_wah_compressed_in_place_and(void)
void test_wah_compressed_in_place_and(void)
{
    /*
     * X
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 1110100000000010101010000000000
     * 0000010000000000000000010101010
     * 00000
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * WAH
     * 536870912
     * 3221225474
     * 1946244096
     * 33554602
     * 0
     *
     * Y
     * |--31-------------------------|
     * 0100000000000000000000000000000
     * 1111111111111111111111111111111
     * 1111111111111111111111111111111
     * 0000000000000000000000000000000
     * 0000000000000000000000000000000
     * 01011
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     * WAH
     * 536870912
     * 3221225474
     * 2147483650
     * 738197504
     *
     *
     * Result:
     * 00100000000000000000000000000000
     * 11000000000000000000000000000010
     * 01110100000000010101010000000000
     * 00000010000000000000000010101010
     * 00101100000000000000000000000000
     */

     unsigned int A[5] = 
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")
        };
 

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };

    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_and(&r_X, &r_Y, &Z);

    unsigned int i;

    // init the first value to be a fill of zeros as long as the max size is
    unsigned int R[6] = {0xc0000006,0,0,0,0,0};


    /*
     * - 00100000000000000000000000000000  - 11000000000000000000000000000110
     *   11000000000000000000000000000010    0
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    0
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000 
     * - 11000000000000000000000000000010  - 11000000000000000000000000000101
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    0
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000 
     *   11000000000000000000000000000010    11000000000000000000000000000010
     * - 10000000000000000000000000000010    0
     *   00101100000000000000000000000000  - 11000000000000000000000000000011
     *                                       0
     *                                       0
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000 
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   10000000000000000000000000000010    0
     * - 00101100000000000000000000000000    10000000000000000000000000000010
     *                                       0
     *                                     - 01111111111111111111111111111111
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000 
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   10000000000000000000000000000010    0
     *   00101100000000000000000000000000    10000000000000000000000000000010
     *                                       0
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000 -> 536870912
     *   11000000000000000000000000000010 -> 3221225474
     *   0                                -> 0
     *   10000000000000000000000000000010 -> 2147483650
     *   0                                -> 0
     *   00101100000000000000000000000000 -> 738197504
     *
     */
    
    unsigned int A0[6] = {536870912,
                          3221225474,
                          0,
                          2147483650,
                          0,
                          738197504};

    unsigned int r = wah_compressed_in_place_and(R, 6, w_Y, wah_size_Y);

    TEST_ASSERT_EQUAL(6, r);
        
    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(A0[i], R[i]);

    r = wah_compressed_in_place_and(R, 6, w_X, wah_size_X);

    /*
     * - 00100000000000000000000000000000  - 00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     * 
     *   00100000000000000000000000000000    00100000000000000000000000000000
     * - 11000000000000000000000000000010  - 11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     * - 01110100000000010101010000000000    0
     *   00000010000000000000000010101010  - 10000000000000000000000000000010
     *   00000000000000000000000000000000    0
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     * - 00000010000000000000000010101010    00000000000000000000000000000000
     *   00000000000000000000000000000000  - 00000000000000000000000000000000
     *                                       00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    00000000000000000000000000000000
     * - 00000000000000000000000000000000    00000000000000000000000000000000
     *                                     - 00101100000000000000000000000000
     *
     *   00100000000000000000000000000000    00100000000000000000000000000000
     *   11000000000000000000000000000010    11000000000000000000000000000010
     *   01110100000000010101010000000000    0
     *   00000010000000000000000010101010    00000000000000000000000000000000
     *   00000000000000000000000000000000    00000000000000000000000000000000
     *                                       00000000000000000000000000000000
     *
     *   00100000000000000000000000000000 -> 536870912
     *   11000000000000000000000000000010 -> 3221225474
     *   0                                -> 0
     *   00000000000000000000000000000000 -> 0
     *   00000010000000000000000010101010 -> 0
     *   00000000000000000000000000000000 -> 0
     *
     */

    unsigned int A1[6] = {536870912,
                          3221225474,
                          0,
                          0,
                          0,
                          0};

    TEST_ASSERT_EQUAL(6, r);
        
    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(A1[i], R[i]);


    unsigned int *ints_o, *ints_ip;
    unsigned int ints_ip_len = compressed_in_place_wah_to_ints(R,r,&ints_ip);
    unsigned int ints_o_len = wah_to_ints(Z,Z_len,&ints_o);

    TEST_ASSERT_EQUAL(ints_o_len, ints_ip_len);
    for (i = 0; i < ints_o_len; ++i)
        TEST_ASSERT_EQUAL(ints_o[i], ints_ip[i]);
}
//}}}

//{{{ void test_wah_in_place_and(void)
void test_wah_in_place_and(void)
{
    /*
     * X
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111111
     * 01000000000101010100000000000000
     * 01000000000000000001010101000000
     * Y
     * |--32--------------------------|
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000001011
     *
     * Result:
     * 01000000000000000000000000000001
     * 11111111111111111111111111111111
     * 11111111111111111111111111111000
     * 00000000000000000000000000000000
     * 00000000000000000000000000000000
     */

    unsigned int A[5] = 
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")
        };

    unsigned int X[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("01000000000101010100000000000000"),
          bin_char_to_int("01000000000000000001010101000000")
        };

    unsigned int Y[5] =
        { bin_char_to_int("01000000000000000000000000000001"),
          bin_char_to_int("11111111111111111111111111111111"),
          bin_char_to_int("11111111111111111111111111111000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000001011")
        };


    unsigned int *w_X;
    int wah_size_X = ints_to_wah(X,5,160,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,160,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_and(&r_X, &r_Y, &Z);

    unsigned int *ints;
    unsigned int ints_len = wah_to_ints(Z, Z_len, &ints);

    unsigned int i;
    for (i = 0; i < 5; ++i) 
        TEST_ASSERT_EQUAL(A[i], ints[i]);


    
    unsigned int R[6] = {0x7fffffff,
                         0x7fffffff,
                         0x7fffffff,
                         0x7fffffff,
                         0x7fffffff,
                         0x7fffffff};

    unsigned int r = wah_in_place_and(R, 6, w_Y, wah_size_Y);
    r = wah_in_place_and(R, 6, w_X, wah_size_X);

    unsigned int *ints_ip;
    unsigned int ints_ip_len = wah_to_ints(R,r,&ints_ip);

    TEST_ASSERT_EQUAL(ints_len, ints_ip_len);
    for (i = 0; i < ints_len; ++i)
        TEST_ASSERT_EQUAL(ints[i], ints_ip[i]);

}
//}}}

//{{{ void test_gt_records_plt_ubin_in_place_wahbm(void)
void test_gt_records_plt_ubin_in_place_wahbm(void)
{
    char *plt_file_name="../data/10.1e4.ind.txt";
    char *ubin_file_name="../data/10.1e4.ind.ubin";
    char *wah_file_name="../data/10.1e4.ind.wahbm";

    struct plt_file pf = init_plt_file(plt_file_name);
    struct ubin_file uf = init_ubin_file(ubin_file_name);
    struct wah_file wf = init_wahbm_file(wah_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *pf_R;
    unsigned int len_pf_R = gt_records_plt(pf, test_records, 4, 0, &pf_R);

    unsigned int *uf_R;
    unsigned int len_uf_R = gt_records_ubin(uf, test_records, 4, 0, &uf_R);

    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_in_place_wahbm(wf,
                                                      test_records,
                                                      4,
                                                      0,
                                                      &wf_R);

    unsigned int *ints;
    unsigned int ints_size = wah_to_ints(wf_R,len_wf_R,&ints);

                
    /*
     * 1000001110111010100101011110111111010110100
     * 0000000220222020000202020000020022000000010
     * 0000000221222021000202010000020022000000000
     * 1011010110222021000202020010020022000000100
     *
     * 0000000110111010000101010000010011000000000
     *
     * 00000001101110100001010100000100 -> 28972292
     * 11000000000                      -> 1536
     */

    unsigned int A[2] = {28972292,1536};
    unsigned int shift[2] = {0,21};
    unsigned int i;
    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
}
//}}}

//{{{ void test_more_wah_compressed_in_place_or(void)
void test_more_wah_compressed_in_place_or(void)
{
    /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.1 same length
     *  1.2 X is longer than Y
     *      1.2.1 X has a single word left
     *      1.2.2 X has a fill left
     *      1.2.2 X has multiple fills left
     *  1.3 Y is longer than X
     *   1.3.1 X is literals
     *   1.3.2 X is fills
     */

    /*
     * 1. Both X and Y are fills
     *  1.1 same length
     */
    unsigned int i;
    unsigned int T1A[2] =
        { bin_char_to_int("11000000000000000000000000000011"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int T1B[2] =
        { bin_char_to_int("11000000000000000000000000000011"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int A1[6] =
        { bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")};


    unsigned int R1[6] = {0x80000006,0,0,0,0,0};
    unsigned int r = wah_compressed_in_place_or(R1, 6, T1A, 2);
    r = wah_compressed_in_place_or(R1, 6, T1B, 2);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A1[i], R1[i]);
    /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.2 X is longer than Y
     *      1.2.1 X has a single word left
     */

    unsigned int T2A[2] =
        { bin_char_to_int("10000000000000000000000000000011"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int T2B[3] =
        { bin_char_to_int("11000000000000000000000000000010"),
         bin_char_to_int("01010101010101010101010101010101"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int A2[6] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")};


    unsigned int R2[6] = {0x80000006,0,0,0,0,0};
    r = wah_compressed_in_place_or(R2, 6, T2A, 2);
    r = wah_compressed_in_place_or(R2, 6, T2B, 3);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A2[i], R2[i]);

     /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.2 X is longer than Y
     *      1.2.2 X has a fill left
     */

    unsigned int T3A[2] =
        { bin_char_to_int("10000000000000000000000000000100"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int T3B[4] =
        { bin_char_to_int("11000000000000000000000000000010"),
         bin_char_to_int("01010101010101010101010101010101"),
         bin_char_to_int("01010101010101010101010101010101"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int A3[7] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")};


    unsigned int R3[7] = {0x80000007,0, 0,0,0,0,0};
    r = wah_compressed_in_place_or(R3, 7, T3A, 2);
    r = wah_compressed_in_place_or(R3, 7, T3B, 4);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A3[i], R3[i]);

     /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.2 X is longer than Y
     *      1.2.3 X has multiple fills left
     */

    unsigned int T4A[1] =
        { bin_char_to_int("10000000000000000000000000000110") };

    unsigned int T4B[3] =
        { bin_char_to_int("11000000000000000000000000000010"),
         bin_char_to_int("01010101010101010101010101010101"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int A4[6] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")};


    unsigned int R4[6] = {0x80000006,0,0,0,0,0};
    r = wah_compressed_in_place_or(R4, 6, T4A, 1);
    r = wah_compressed_in_place_or(R4, 6, T4B, 3);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A4[i], R4[i]);

     /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.3 Y is longer than X
     *   1.3.1 X is literals
     */

    unsigned int T5A[5] =
        { bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("01111111111111111111111111111111"),
          bin_char_to_int("00000000000000000000000000000001"),
          bin_char_to_int("01110001110001110001110001110001"),
          bin_char_to_int("01110001110001110001110001110001")
        };

    unsigned int T5B[2] =
        { bin_char_to_int("10000000000000000000000000000100"),
          bin_char_to_int("00001110001110001110001110001110")};

    unsigned int A5[5] =
        { bin_char_to_int("01010101010101010101010101010101"),
          bin_char_to_int("01111111111111111111111111111111"),
          bin_char_to_int("00000000000000000000000000000001"),
          bin_char_to_int("01110001110001110001110001110001"),
          bin_char_to_int("01111111111111111111111111111111")
        };

    unsigned int R5[6] = {0x80000005,0,0,0,0};
    r = wah_compressed_in_place_or(R5, 5, T5A, 5);
    r = wah_compressed_in_place_or(R5, 5, T5B, 2);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A5[i], R5[i]);

     /*
     * Cases:
     * 1. Both X and Y are fills
     *  1.3 Y is longer than X
     *   1.3.2 X is fills
     */

    unsigned int T6A[3] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("10000000000000000000000000000010"),
          bin_char_to_int("11000000000000000000000000000010"),
        };

    unsigned int T6B[2] =
        { bin_char_to_int("10000000000000000000000000000100"),
          bin_char_to_int("11000000000000000000000000000010")};

    unsigned int A6[6] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("10000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000")
        };

    unsigned int R6[6] = {0x80000006,0,0,0,0,0};
    r = wah_compressed_in_place_or(R6, 6, T6A, 3);
    r = wah_compressed_in_place_or(R6, 6, T6B, 2);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A6[i], R6[i]);

    // X is a fill of ones with a litteral left over
    //
    unsigned int T7A[2] =
        { bin_char_to_int("11000000000000000000000000000011"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int T7B[3] =
        { bin_char_to_int("11000000000000000000000000000010"),
         bin_char_to_int("01010101010101010101010101010101"),
         bin_char_to_int("11000000000000000000000000000011")};

    unsigned int A7[6] =
        { bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("01111111111111111111111111111111"),
          bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000")};


    unsigned int R7[6] = {0x80000006,0,0,0,0,0};
    r = wah_compressed_in_place_or(R7, 6, T7A, 2);
    r = wah_compressed_in_place_or(R7, 6, T7B, 3);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A7[i], R7[i]);


    // Y is a fill of 1s, and the last value 

    unsigned int T8A[3] =
        { bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000011"),
          bin_char_to_int("10000000000000000000000000000010")
        };

    unsigned int T8B[2] =
        { bin_char_to_int("11000000000000000000000000000100"),
          bin_char_to_int("11000000000000000000000000000010")};

    unsigned int A8[6] =
        { bin_char_to_int("11000000000000000000000000000011"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("00000000000000000000000000000000"),
          bin_char_to_int("01111111111111111111111111111111"),
          bin_char_to_int("11000000000000000000000000000010"),
          bin_char_to_int("00000000000000000000000000000000")
        };

    unsigned int R8[6] = {0x80000006,0,0,0,0,0};
    r = wah_compressed_in_place_or(R8, 6, T8A, 3);
    r = wah_compressed_in_place_or(R8, 6, T8B, 2);

    for (i = 0; i < r; ++i)
        TEST_ASSERT_EQUAL(A8[i], R8[i]);


}
//}}}

//{{{void test_plt_gt_count(void)
void test_plt_gt_count(void)
{
    /*
     * Input:
     * 1000001110111010100101011110111111010110100
     * 0000000220222020000202020000020022000000010
     * 0000000221222021000202010000020022000000000
     * 1011010110222021000202020010020022000000100
     *
     */

    unsigned int A[43] = {
        2,0,1,1,0,1,1,4,4,1,4,4,4,0,4,2,1,
        0,0,4,0,4,0,4,1,1,2,0,1,4,1,1,4,4,
        0,1,0,1,1,0,2,1,0};

    char *plt_file_name="../data/10.1e4.ind.txt";
    struct plt_file pf = init_plt_file(plt_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *R;
    unsigned int R_len = count_range_records_plt(pf,
                                                 test_records,
                                                 4,
                                                 0,
                                                 3,
                                                 &R);
    unsigned int i;
    for (i = 0; i < R_len; ++i)
        TEST_ASSERT_EQUAL(A[i], R[i]);

}
//}}}

//{{{void test_plt_gt_count(void)
void test_ubin_gt_count(void)
{
#if 0
    /*
     * Input:
     * 1000001110111010100101011110111111010110100
     * 0000000220222020000202020000020022000000010
     * 0000000221222021000202010000020022000000000
     * 1011010110222021000202020010020022000000100
     *
     */

    unsigned int A[43] = {
        2,0,1,1,0,1,1,4,4,1,4,4,4,0,4,2,1,
        0,0,4,0,4,0,4,1,1,2,0,1,4,1,1,4,4,
        0,1,0,1,1,0,2,1,0};

    char *ubin_file_name="../data/10.1e4.ind.ubin";

    struct ubin_file uf = init_ubin_file(ubin_file_name);

    unsigned int test_records[4] = {1,2,3,4};

    unsigned int *R;
    unsigned int R_len = count_range_records_ubin(uf,
                                                 test_records,
                                                 4,
                                                 0,
                                                 3,
                                                 &R);
    unsigned int i;
    for (i = 0; i < R_len; ++i)
        TEST_ASSERT_EQUAL(A[i], R[i]);
#endif
}
//}}}

//{{{void test add_wahmb(void)
void test_add_wahmb(void)
{
    unsigned int wah1[4] = {
        bin_char_to_int("11000000000000000000000000000010"), //2
        bin_char_to_int("01010101010101010101010101010101"), //1
        bin_char_to_int("10000000000000000000000000000011"), //3
        bin_char_to_int("01010101010101010101010101010101") //1
    };

    unsigned int i, *R, r_size;
    r_size = (2+1+3+1)*31;
    R = (unsigned int*) calloc(r_size, sizeof(unsigned int));

    r_size = add_wahbm(R, r_size, wah1, 4);
    r_size = add_wahbm(R, r_size, wah1, 4);
    r_size = add_wahbm(R, r_size, wah1, 4);

    unsigned int sum = 0;
    for (i = 0; i < r_size; ++i)
        sum += R[i];

    TEST_ASSERT_EQUAL((2*31+16*2)*3,sum);

}
//}}}

//{{{void test add_compressed_in_place_wahmb(void)
void test_add_compressed_in_place_wahmb(void)
{
    unsigned int wah1[7] = {
        bin_char_to_int("11000000000000000000000000000010"), //2
        bin_char_to_int("00000000000000000000000000000000"), //blank
        bin_char_to_int("01010101010101010101010101010101"), //1
        bin_char_to_int("10000000000000000000000000000011"), //3
        bin_char_to_int("00000000000000000000000000000000"), //blank
        bin_char_to_int("00000000000000000000000000000000"), //blank
        bin_char_to_int("01010101010101010101010101010101") //1
    };

    unsigned int i, *R, r_size;
    r_size = (7)*31;
    R = (unsigned int*) calloc(r_size, sizeof(unsigned int));

    r_size = add_compressed_in_place_wahbm(R, r_size, wah1, 7);
    r_size = add_compressed_in_place_wahbm(R, r_size, wah1, 7);
    r_size = add_compressed_in_place_wahbm(R, r_size, wah1, 7);

    unsigned int sum = 0;
    for (i = 0; i < r_size; ++i) 
        sum += R[i];

    TEST_ASSERT_EQUAL((2*31+16*2)*3,sum);

}
//}}}

//{{{void test_invert_plt_ubin(void)
void test_invert_plt_ubin(void)
{
    char *l[10] = {
        "2 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 2 1 0 0 1 0 1 0 1 1 0 1 1 "
            "1 1 1 1 0 1 1 1 1 1 1 0 1 0",
        "1 0 0 0 0 0 1 1 1 0 1 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 1 0 1 "
            "1 1 1 1 1 0 1 0 1 1 0 1 0 0",
        "0 0 0 0 0 0 0 2 2 0 2 2 2 0 2 0 0 0 0 2 0 2 0 2 0 0 0 0 0 "
            "2 0 0 2 2 0 0 0 0 0 0 0 1 0",
        "0 0 0 0 0 0 0 2 2 1 2 2 2 0 2 1 0 0 0 2 0 2 0 1 0 0 0 0 0 "
            "2 0 0 2 2 0 0 0 0 0 0 0 0 0",
        "1 0 1 1 0 1 0 1 1 0 2 2 2 0 2 1 0 0 0 2 0 2 0 2 0 0 1 0 0 "
            "2 0 0 2 2 0 0 0 0 0 0 1 0 0",
        "0 0 0 0 0 0 0 2 2 0 2 2 2 0 2 0 0 0 1 2 0 2 0 2 0 0 0 0 0 "
            "2 0 0 2 2 0 0 0 0 0 0 0 0 0",
        "1 0 0 0 2 0 2 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 "
            "2 0 0 2 2 0 0 0 0 0 0 0 0 0",
        "0 0 0 0 0 0 0 2 2 2 1 1 1 0 1 1 1 0 0 1 0 1 0 0 1 1 0 0 1 "
            "1 1 0 1 0 1 1 1 1 1 0 0 1 1",
        "0 1 0 0 0 0 0 2 2 0 1 1 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 0 "
            "2 0 0 2 1 1 0 1 0 0 1 1 1 1",
        "1 0 1 0 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 0 0 0 0 "
            "2 0 0 2 1 1 0 1 0 0 0 0 0 0",
    };

    unsigned int A[43][10] = {
            {2,1,0,0,1,0,1,0,0,1},
            {0,0,0,0,0,0,0,0,1,0},
            {1,0,0,0,1,0,0,0,0,1},
            {1,0,0,0,1,0,0,0,0,0},
            {0,0,0,0,0,0,2,0,0,0},
            {1,0,0,0,1,0,0,0,0,1},
            {1,1,0,0,0,0,2,0,0,0},
            {0,1,2,2,1,2,0,2,2,1},
            {0,1,2,2,1,2,0,2,2,1},
            {0,0,0,1,0,0,0,2,0,0},
            {0,1,2,2,2,2,0,1,1,1},
            {0,1,2,2,2,2,0,1,1,1},
            {0,1,2,2,2,2,0,1,1,1},
            {1,0,0,0,0,0,0,0,0,1},
            {0,1,2,2,2,2,0,1,1,1},
            {0,0,0,1,1,0,0,1,0,0},
            {2,1,0,0,0,0,2,1,1,1},
            {1,0,0,0,0,0,0,0,0,1},
            {0,0,0,0,0,1,0,0,0,0},
            {0,1,2,2,2,2,0,1,1,1},
            {1,0,0,0,0,0,0,0,0,1},
            {0,1,2,2,2,2,0,1,1,1},
            {1,0,0,0,0,0,0,0,0,1},
            {0,1,2,1,2,2,0,0,1,1},
            {1,1,0,0,0,0,0,1,0,0},
            {1,1,0,0,0,0,0,1,0,0},
            {0,1,0,0,1,0,0,0,1,0},
            {1,0,0,0,0,0,0,0,0,0},
            {1,1,0,0,0,0,0,1,0,0},
            {1,1,2,2,2,2,2,1,2,2},
            {1,1,0,0,0,0,0,1,0,0},
            {1,1,0,0,0,0,0,0,0,0},
            {1,1,2,2,2,2,2,1,2,2},
            {0,1,2,2,2,2,2,0,1,1},
            {1,0,0,0,0,0,0,1,1,1},
            {1,1,0,0,0,0,0,1,0,0},
            {1,0,0,0,0,0,0,1,1,1},
            {1,1,0,0,0,0,0,1,0,0},
            {1,1,0,0,0,0,0,1,0,0},
            {1,0,0,0,0,0,0,0,1,0},
            {0,1,0,0,1,0,0,0,1,0},
            {1,0,1,0,0,0,0,1,1,0},
            {0,0,0,0,0,0,0,1,1,0}
    };

    unsigned int o_num_fields, o_num_records;
    o_num_fields = 43;
    o_num_records = 10;

    unsigned int n_num_fields, n_num_records;
    unsigned int two_bit_i = 0;
    unsigned int **ubin = NULL;

    unsigned int i,j;

    for (i = 0; i < o_num_records; ++i)
        two_bit_i = invert_plt_to_ubin(l[i],
                                    o_num_fields,
                                    o_num_records,
                                    &n_num_fields,
                                    &n_num_records,
                                    two_bit_i,
                                    &ubin);

    TEST_ASSERT_EQUAL(o_num_fields, n_num_records);
    TEST_ASSERT_EQUAL(o_num_records, n_num_fields);

    for (i = 0; i < n_num_records; ++i) {
        int *r = unpack_2_bit_ints(ubin[i][0]);
        for (j = 0; j < n_num_fields; ++j) 
            TEST_ASSERT_EQUAL(A[i][j],r[j]);
        free(r);
    }
}
//}}}

//{{{ void test_convert_file_by_name_invert_plt_to_ubin(void)
void test_convert_file_by_name_invert_plt_to_ubin(void)
{

    char *plt_file = "../data/10.1e4.ind.txt";
    char *i_ubin_file = "tmp.i.ubin";

    convert_file_by_name_invert_plt_to_ubin(plt_file, i_ubin_file);
}
//}}}

//{{{void test_get_wah_bitmaps_in_place(void)
void test_get_wah_bitmaps_in_place(void)
{
    char *wahbm_file_name = "../data/10.1e4.ind.wahbm";
    struct wah_file wf = init_wahbm_file(wahbm_file_name);

    unsigned int bm_size_0, bm_size_1, bm_size_2, bm_size_3;

    unsigned int *bm_0,*bm_1, *bm_2, *bm_3;

    unsigned int max_wah_size = (wf.num_fields + 31 - 1)/ 31;

    bm_0 = (unsigned int *) malloc(sizeof(unsigned int)*max_wah_size);
    bm_1 = (unsigned int *) malloc(sizeof(unsigned int)*max_wah_size);
    bm_2 = (unsigned int *) malloc(sizeof(unsigned int)*max_wah_size);
    bm_3 = (unsigned int *) malloc(sizeof(unsigned int)*max_wah_size);

    bm_size_0 = get_wah_bitmap_in_place(wf, 1, 0, &bm_0);
    bm_size_1 = get_wah_bitmap_in_place(wf, 1, 1, &bm_1);
    bm_size_2 = get_wah_bitmap_in_place(wf, 1, 2, &bm_2);
    bm_size_3 = get_wah_bitmap_in_place(wf, 1, 3, &bm_3);

    unsigned int bms_size;
    unsigned int bms_sizes[4];
    unsigned int *bms;
    bms = (unsigned int *) malloc(sizeof(unsigned int)*max_wah_size*4);
    bms_size = get_wah_bitmaps_in_place(wf, 1, &bms, bms_sizes);

    TEST_ASSERT_EQUAL(bm_size_0+ bm_size_1+ bm_size_2+ bm_size_3, bms_size);

    TEST_ASSERT_EQUAL(bm_size_0,bms_sizes[0]);
    TEST_ASSERT_EQUAL(bm_size_1,bms_sizes[1]);
    TEST_ASSERT_EQUAL(bm_size_2,bms_sizes[2]);
    TEST_ASSERT_EQUAL(bm_size_3,bms_sizes[3]);


    unsigned int *bms_i = bms;
    unsigned int i;

    for (i = 0; i < bm_size_0; ++i)
        TEST_ASSERT_EQUAL(bm_0[i], bms_i[i]);

    bms_i = bms + bms_sizes[0];

    for (i = 0; i < bm_size_1; ++i)
        TEST_ASSERT_EQUAL(bm_1[i], bms_i[i]);

    bms_i += bms_sizes[1];

    for (i = 0; i < bm_size_2; ++i)
        TEST_ASSERT_EQUAL(bm_2[i], bms_i[i]);

    bms_i += bms_sizes[2];

    for (i = 0; i < bm_size_3; ++i)
        TEST_ASSERT_EQUAL(bm_3[i], bms_i[i]);

}
//}}}

//{{{ void test_gt_records_fields(void)
void test_gt_records_fields(void)
{
    char *plt_ind_file_name="../data/10.1e4.ind.txt";
    char *plt_var_file_name="../data/10.1e4.var.txt";

    struct plt_file pf_i = init_plt_file(plt_ind_file_name);
    struct plt_file pf_v = init_plt_file(plt_var_file_name);

    unsigned int test_records[4] = {1,2,3,4};
    unsigned int test_fields[4] = {1,2,3,4};

    unsigned int *pf_R_i;
    unsigned int len_pf_R_i = gt_records_plt(pf_i, test_records, 4, 0, &pf_R_i);

    unsigned int *pf_R_v;
    unsigned int len_pf_R_v = gt_fields_plt(pf_v, test_fields, 4, 0, &pf_R_v);

    /*
    printf("%u %u \n", len_pf_R_i, len_pf_R_v);
    unsigned int i;
    for (i = 0; i < 2; ++i)
        printf("%u %u \n", pf_R_i[i], pf_R_v[i]);

    unsigned int A[2] = {28972292,1536};
    unsigned int shift[2] = {0,21};
    unsigned int i;
    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , pf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , uf_R[i] >> shift[i]);

    for (i = 0; i < 2; ++i)
        TEST_ASSERT_EQUAL(A[i] , ints[i] >> shift[i]);
    */
}
//}}}

//{{{void test_avx_gt_count_records_in_place_wahbm(void)
void test_avx_gt_count_records_in_place_wahbm(void)
{
    char *in = "../data/10.1e4.ind.wahbm";
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R;
    unsigned int num_records = 5;
    unsigned int R[5] = {1,2,3,4,5};
    unsigned int query_value = 0;

    len_wf_R = gt_count_records_in_place_wahbm(wf,
                                               R,
                                               num_records,
                                               query_value,
                                               &wf_R);

    unsigned int avx_len_wf_R;
    unsigned int *avx_wf_R;
    avx_len_wf_R = avx_gt_count_records_in_place_wahbm(wf,
                                                   R,
                                                   num_records,
                                                   query_value,
                                                   &avx_wf_R);

    TEST_ASSERT_EQUAL(len_wf_R,avx_len_wf_R);

    unsigned int i;

    for (i = 0; i < len_wf_R; ++i) 
        TEST_ASSERT_EQUAL(wf_R[i], avx_wf_R[i]);

}
//}}}

//{{{void test_avx_add_wahbm(void)
void test_avx_add_wahbm(void)
{
    uint32_t W1[5] = {
        bin_char_to_int("01111000000000000000000000001111"), // 1
        bin_char_to_int("10000000000000000000000000000010"), // 2 
        bin_char_to_int("01110010101010101010101010100111"), // 1
        bin_char_to_int("01100010101010101010101010100011"), // 1
        bin_char_to_int("11000000000000000000000000000010")  // 2
    };

    uint32_t A1[7*31] = {
    3,3,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,
    2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    2,2,2,1,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1,2,2,2,
    2,2,1,1,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1,1,2,2,
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

    unsigned int *R1;
    int r = posix_memalign((void **)&R1, 32, 7 *31* sizeof(unsigned int));

    unsigned int R1_len = avx_add_wahbm(R1,
                                        7*31,
                                        W1,
                                        5);

    unsigned int i;

    uint32_t W2[1] = {
        bin_char_to_int("11000000000000000000000000000111") //7
    };

    R1_len = avx_add_wahbm(R1,
                           7*31,
                           W2,
                           1);

    uint32_t W3[3] = {
        bin_char_to_int("01111110000000000000000000111111"), //1
        bin_char_to_int("01110000000000000000000000000111"), //1
        bin_char_to_int("10000000000000000000000000000101") //5
    };

    R1_len = avx_add_wahbm(R1,
                           7*31,
                           W3,
                           3);

    /*
    for (i = 0; i < R1_len; ++i) {
        if ( (i>0) && (i%31==0))
            printf("\n");
        printf("%u", R1[i]);
    }
    printf("\n");
    */


    for (i = 0; i < R1_len; ++i) 
        TEST_ASSERT_EQUAL(A1[i], R1[i]);
}
//}}}

//{{{void test_avx_add_wahbm(void)
void test_avx_add_wahbm_2(void)
{
    uint32_t W1[3] = { 0 , 0 , 0};
    uint32_t W2[3] = { 0 , 0 , 0};
    uint32_t W3[3] = { 8192 , 0 , 2101248};
    uint32_t W4[3] = { 0 , 0 , 0};
    uint32_t W5[3] = { 0 , 512 , 0};

    unsigned int *R1;
    int r = posix_memalign((void **)&R1, 32, 3 *31* sizeof(unsigned int));

    unsigned int i;
    unsigned int R1_len = avx_add_wahbm(R1, 3*31, W1, 3);
    R1_len = avx_add_wahbm(R1, 3*31, W2, 3);
    R1_len = avx_add_wahbm(R1, 3*31, W3, 3);
    R1_len = avx_add_wahbm(R1, 3*31, W4, 3);
    R1_len = avx_add_wahbm(R1, 3*31, W5, 3);

    uint32_t A[93] = { 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0
    };


    for (i = 0; i < R1_len; ++i) 
        TEST_ASSERT_EQUAL(A[i], R1[i]);
}
//}}}

//{{{void test_pq(void)
void test_pq(void)
{

    pri_queue q = priq_new(0);

    priority p1, p2, p3, p4, p5;

    p1.sum = 3;
    p1.len = 5;
    int p1_d = 1;

    p2.sum = 4;
    p2.len = 2;
    int p2_d = 2;

    p3.sum = 4;
    p3.len = 1;
    int p3_d = 3;

    p4.sum = 4;
    p4.len = 3;
    int p4_d = 4;

    p5.sum = 1;
    p5.len = 1;
    int p5_d = 5;

    priq_push(q, &p1_d, p1);
    priq_push(q, &p2_d, p2);
    priq_push(q, &p3_d, p3);
    priq_push(q, &p4_d, p4);
    priq_push(q, &p5_d, p5);

    priority p;

    uint32_t A[5] = {5, 1, 4, 2, 3};
    uint32_t a_i = 0;

    while ( priq_top(q, &p) != NULL ) {
        int *d = priq_pop(q, &p);
        TEST_ASSERT_EQUAL(A[a_i], *d);
        a_i += 1;
    }
}
//}}}

//{{{void test_push_bcf_gt_md(void)
void test_push_bcf_gt_md(void)
{
    uint32_t num_inds = 10;
    uint32_t num_vars = 43;
    
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    char *gt_of_name = ".gt.tmp.packed";
    char *md_of_name = ".md.tmp.packed";

    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);
    pri_queue q = priq_new(0);

    uint32_t *md_index = (uint32_t *) malloc(num_vars * sizeof(uint32_t));

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);


    // some select entries from the meta data
    char *md_A_0 = "1\t1\tV1\tA\tT";
    char *md_A_1 = "1\t2\tV2\tA\tT";
    char *md_A_30 = "1\t31\tV31\tA\tT";
    char *md_A_42 = "1\t43\tV43\tA\tT";

    FILE *md_of = fopen(md_of_name, "r");

    // Test 0
    uint32_t start = 0;
    uint32_t len = md_index[0] - start;
    fseek(md_of, start*sizeof(char), SEEK_SET);
    char buf_0[len+1];
    fread(buf_0, sizeof(char), len, md_of);
    buf_0[len] = '\0';
    TEST_ASSERT_EQUAL(0, strcmp(buf_0, md_A_0));

    // Test 1
    start = md_index[1-1];
    len = md_index[1] - start;
    fseek(md_of, start*sizeof(char), SEEK_SET);
    char buf_1[len+1];
    fread(buf_1, sizeof(char), len, md_of);
    buf_1[len] = '\0';
    TEST_ASSERT_EQUAL(0, strcmp(buf_1, md_A_1));

    // Test 30
    start = md_index[30-1];
    len = md_index[30] - start;
    fseek(md_of, start*sizeof(char), SEEK_SET);
    char buf_30[len+1];
    fread(buf_30, sizeof(char), len, md_of);
    buf_30[len] = '\0';
    TEST_ASSERT_EQUAL(0, strcmp(buf_30, md_A_30));

    // Test 42
    start = md_index[42-1];
    len = md_index[42] - start;
    fseek(md_of, start*sizeof(char), SEEK_SET);
    char buf_42[len+1];
    fread(buf_42, sizeof(char), len, md_of);
    buf_42[len] = '\0';
    TEST_ASSERT_EQUAL(0, strcmp(buf_42, md_A_42));

    fclose(md_of);
    remove(md_of_name);

    /*
     * Test some select lines from the genotypes:
     * 0:  2 1 0 0 1 0 1 0 0 1 -> 2420379648
     * 1:  0 0 0 0 0 0 0 0 1 0 -> 16384
     * 30: 1 1 0 0 0 0 0 1 0 0 -> 1342242816
     * 42: 0 0 0 0 0 0 0 1 1 0 -> 81920
     */
    uint32_t gt_A[4] = {2420379648, 16384, 1342242816, 81920};
    uint32_t gt_T[4] = {0,1,30,42};

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16); // should be 1

    FILE *gt_of = fopen(gt_of_name, "rb");

    uint32_t i_buff;

    uint32_t i;
    for (i = 0; i < 4; ++i) {
        fseek(gt_of, gt_T[i]*num_ind_ints*sizeof(uint32_t), SEEK_SET);
        fread(&i_buff, sizeof(uint32_t), num_ind_ints, gt_of);
        TEST_ASSERT_EQUAL(gt_A[i], i_buff);
    }
    fclose(gt_of);
    remove(gt_of_name);

    // Test the priority q order
    priority p;
    i = 0;
    uint32_t last_sum = 0, last_len = 10;
    while ( priq_top(q, &p) != NULL ) {
        int *d = priq_pop(q, &p);
        if (p.sum == last_sum)
            TEST_ASSERT_EQUAL(1, last_len >= p.len);
        else
            TEST_ASSERT_EQUAL(1, last_sum <= p.sum);
        last_sum = p.sum;
        last_len = p.len;
    }

    priq_free(q);
    close_bcf_file(&bcf_f); 
}
//}}}

//{{{void test_sort_gt_md(void)
void test_sort_gt_md(void)
{
    uint32_t num_inds = 10;
    uint32_t num_vars = 43;
    
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    char *gt_of_name = ".gt.tmp.packed";
    char *s_gt_of_name = ".s.gt.tmp.packed";
    char *md_of_name = ".md.tmp.packed";
    char *bim_name = "md.bim";
    char *vid_name = "vid";

    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);
    pri_queue q_t = priq_new(0);
    pri_queue q = priq_new(0);

    uint32_t *md_index = (uint32_t *) malloc(num_vars * sizeof(uint32_t));

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    sort_gt_md(&q,
               md_index,
               num_inds,
               num_vars,
               gt_of_name,
               s_gt_of_name,
               md_of_name,
               bim_name,
               vid_name);

    // close the bcf and open/process it again so we can get the 
    // prioirity q
    close_bcf_file(&bcf_f); 
    bcf_f = init_bcf_file(bcf_file_name);
    push_bcf_gt_md(&q_t,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    // Push the zero-based ids into an arry
    uint32_t P[43];
    uint32_t i = 0;
    priority p;
    while ( priq_top(q_t, &p) != NULL ) {
        int *d = priq_pop(q_t, &p);
        P[i] = *d;
        i+=1;
    }

    // Make sure the zero-based ids match the one-based id in the 
    // sorted bim file
    FILE *md_out = fopen(bim_name,"r");
    char line[100];
    i = 0;
    char *pch;
    while(fgets(line, 100, md_out) != NULL) {
        pch = strtok(line, "\t");
        pch = strtok(NULL, "\t");
        TEST_ASSERT_EQUAL(P[i]+1, atoi(pch));
        i+=1;
    }
    fclose(md_out);


    // Grab the sorted file in order and the unsored file in the order
    // defined by the priority queue and make sure they match
    FILE *gt_of = fopen(gt_of_name,"rb");
    FILE *s_gt_of = fopen(s_gt_of_name,"rb");
    uint32_t gt, s_gt;
    for (i = 0; i < num_vars; ++i) {
        fseek(s_gt_of, i*sizeof(uint32_t), SEEK_SET);
        fread(&s_gt, sizeof(uint32_t), 1, s_gt_of);

        fseek(gt_of, P[i]*sizeof(uint32_t), SEEK_SET);
        fread(&gt, sizeof(uint32_t), 1, gt_of);

        TEST_ASSERT_EQUAL(gt, s_gt);
    }
    fclose(gt_of);
    fclose(s_gt_of);

    priq_free(q);
    priq_free(q_t);
    close_bcf_file(&bcf_f); 
    free(md_index);
}
//}}}

//{{{void test_rotate_encode_wahbm(void)
void test_rotate_encode_wahbm(void)
{
    uint32_t num_inds = 10;
    uint32_t num_vars = 43;
    
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    char *gt_of_name = ".gt.tmp.packed";
    char *s_gt_of_name = ".s.gt.tmp.packed";
    char *r_s_gt_of_name = ".r.s.gt.tmp.packed";
    char *md_of_name = ".md.tmp.packed";
    char *bim_name = "md.bim";
    char *vid_name = "vid";

    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);
    pri_queue q = priq_new(0);

    uint32_t *md_index = (uint32_t *) malloc(num_vars * sizeof(uint32_t));

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    sort_gt_md(&q,
               md_index,
               num_inds,
               num_vars,
               gt_of_name,
               s_gt_of_name,
               md_of_name,
               bim_name,
               vid_name);

    rotate_encode_wahbm(num_inds,
                        num_vars,
                        s_gt_of_name,
                        r_s_gt_of_name);

    FILE *s_gt_of = fopen(s_gt_of_name,"rb");
    uint32_t n_r[43];
    uint32_t i;
    for (i = 0; i < num_vars; ++i) {
        fseek(s_gt_of, i*sizeof(uint32_t), SEEK_SET);
        fread(&(n_r[i]), sizeof(uint32_t), 1, s_gt_of);
    }
    fclose(s_gt_of);

    uint32_t num_var_ints = 1 + ((num_vars - 1) / 16);

    FILE *r_s_gt_of = fopen(r_s_gt_of_name,"rb");
    uint32_t *r = (uint32_t *) 
            malloc(num_inds*num_var_ints*sizeof(uint32_t));
    for (i = 0; i < num_vars; ++i) {
        fseek(r_s_gt_of, 2*sizeof(uint32_t), SEEK_SET);
        fread(r, sizeof(uint32_t), num_inds*num_var_ints, r_s_gt_of);
    }
    fclose(r_s_gt_of);


    // test ind 0
    uint32_t v = 0;

    for (i = 0; i < 16; ++i) {
        v += ((n_r[i] >> 30) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[0]);

    v = 0;
    for (; i < 32; ++i) {
        v += ((n_r[i] >> 30) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[1]);

    v = 0;
    for (; i < num_vars; ++i) {
        v += ((n_r[i] >> 30) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[2]);

    // test ind 6 
    v = 0;
    for (i = 0; i < 16; ++i) {
        v += ((n_r[i] >> (30-2*6)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[0+num_var_ints*6]);

    v = 0;
    for (; i < 32; ++i) {
        v += ((n_r[i] >> (30-2*6)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[1+num_var_ints*6]);

    v = 0;
    for (; i < num_vars; ++i) {
        v += ((n_r[i] >> (30-2*6)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[2+num_var_ints*6]);

    // test ind 9
    v = 0;
    for (i = 0; i < 16; ++i) {
        v += ((n_r[i] >> (30-2*9)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[0+num_var_ints*9]);

    v = 0;
    for (; i < 32; ++i) {
        v += ((n_r[i] >> (30-2*9)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[1+num_var_ints*9]);

    v = 0;
    for (; i < num_vars; ++i) {
        v += ((n_r[i] >> (30-2*9)) & 3) << (30 - i*2);
    }
    TEST_ASSERT_EQUAL(v, r[2+num_var_ints*9]);


    priq_free(q);
    close_bcf_file(&bcf_f); 
    free(md_index);
}
//}}}

void test_parse_q(void)
{
    struct gqt_query q;

    char qt[100];
    int r;

    strcpy(qt, "HET");
    r = parse_q(qt, &q);

    TEST_ASSERT_EQUAL(0, r);

    TEST_ASSERT_EQUAL(-1, q.variant_op);
    TEST_ASSERT_EQUAL(-1, q.op_condition);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[0]);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[1]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[2]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[3]);
    TEST_ASSERT_EQUAL(-1, q.condition_value);


    strcpy(qt, "HET HOMO_ALT");
    r = parse_q(qt, &q);

    TEST_ASSERT_EQUAL(0, r);

    TEST_ASSERT_EQUAL(-1, q.variant_op);
    TEST_ASSERT_EQUAL(-1, q.op_condition);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[0]);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[1]);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[2]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[3]);
    TEST_ASSERT_EQUAL(-1, q.condition_value);


    strcpy(qt, "count(HET) >= 3");
    r = parse_q(qt, &q);

    TEST_ASSERT_EQUAL(0, r);

    TEST_ASSERT_EQUAL(p_count, q.variant_op);
    TEST_ASSERT_EQUAL(p_greater_than_equal, q.op_condition);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[0]);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[1]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[2]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[3]);
    TEST_ASSERT_EQUAL(3, q.condition_value);

    strcpy(qt, "pct(HOMO_REF HET) != 0.001");
    r = parse_q(qt, &q);

    TEST_ASSERT_EQUAL(0, r);

    TEST_ASSERT_EQUAL(p_pct, q.variant_op);
    TEST_ASSERT_EQUAL(p_not_equal, q.op_condition);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[0]);
    TEST_ASSERT_EQUAL(1, q.genotype_condition[1]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[2]);
    TEST_ASSERT_EQUAL(0, q.genotype_condition[3]);
    TEST_ASSERT_EQUAL(0.001, q.condition_value);


    // test bad syntax
    strcpy(qt, "count()");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);

    strcpy(qt, "count(HET");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);

    strcpy(qt, "count(HET) >=");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);

    strcpy(qt, "count(HET) 10");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);

    strcpy(qt, "HET == 10");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);

    strcpy(qt, "HET ==");
    r = parse_q(qt, &q);
    TEST_ASSERT_EQUAL(1, r);
}

//{{{  OLD

#if 0
//{{{ void test_rotate_sort_bcf_to_wahbm(void)
//void test_rotate_sort_bcf_to_wahbm(void)
{
    uint32_t num_vars = 43;
    uint32_t num_inds = 10;

    char *hdf5_file_name = "../data/10.1e4.var.h5";
    struct hdf5_file hdf5_f = init_hdf5_file(hdf5_file_name,
                                             num_vars,
                                             num_inds);

    char *bcf_file_name = "../data/10.1e4.var.bcf";
    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);

    pri_queue q = priq_new(0);

    push_bcf_gt_md(&q, &bcf_f, &hdf5_f);

    uint32_t *gt = (uint32_t *) malloc(hdf5_f.num_gt_ints*sizeof(uint32_t));
    int r = read_hdf5_gt(hdf5_f, 3, gt);

    char *md_out;
    read_hdf5_md(hdf5_f, 0, &md_out);

    //1 4   V4  A   T   100 PASS    N=A GT  
    //0|1 0|0 0|0 0|0 0|1 0|0 0|0 0|0 0|0 0|0 
    //1 0 0 0 1 0 0 0 0 0
    //01000000010000000000000000000000 -> 1077936128
    
    TEST_ASSERT_EQUAL(1077936128, gt[0]);
    //TEST_ASSERT_EQUAL(0, strcmp("1\t4\tV4\tA\tT", md_out));


    sort_rotate_gt_md(&q, &hdf5_f, "tmp.bim");

    /*01000001010101010101000001010101 -> 1096110165
     *01010101010000000001010101101000 -> 1430263144
     *00000000000000000001010000000000 -> 5120
     *                      |--------|
     */

    free(gt);
    gt = (uint32_t *) malloc(hdf5_f.num_r_gt_ints*sizeof(uint32_t));
    r = read_hdf5_r_gt(hdf5_f, 0, gt);

    TEST_ASSERT_EQUAL(1096110165, gt[0]);
    TEST_ASSERT_EQUAL(1430263144, gt[1]);
    TEST_ASSERT_EQUAL(5120, gt[2]);


    uint32_t A_0[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    uint32_t A_5[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    uint32_t i;

    for  (i = 0; i < num_inds; ++i){
        r = read_hdf5_r_gt(hdf5_f, i, gt);
        TEST_ASSERT_EQUAL(A_0[i], gt[0] >> 30);
        TEST_ASSERT_EQUAL(A_5[i], (gt[0] >> (30-5*2))&3);
    }

    r = convert_hdf5_ind_ubin_to_ind_wah(hdf5_f, "tmp.hdf5.wah");

    struct wah_file wf = init_wahbm_file("tmp.hdf5.wah");


    TEST_ASSERT_EQUAL(num_inds,wf.num_records);
    TEST_ASSERT_EQUAL(num_vars,wf.num_fields);

    unsigned int test_record, test_bitmap;
    unsigned int *ints, num_ints;
    unsigned int *wah_bms[4], wah_sizes[4];
    unsigned int *wah_ints[4], wah_num_ints[4];
    unsigned int two_bit, bit, ubin_int_i, ubin_bit_i,
                 wah_int_i, wah_bit_i, field_i;

    for (test_record = 0; test_record < 8; ++test_record) {
        field_i = 0;
        //num_ints = get_ubin_record(uf, test_record, &ints);
        r = read_hdf5_r_gt(hdf5_f, test_record, gt);
        num_ints = hdf5_f.num_r_gt_ints;
        for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap) {
            wah_sizes[test_bitmap] = get_wah_bitmap(wf,
                                                    test_record,
                                                    test_bitmap,
                                                    &(wah_bms[test_bitmap]));
            wah_num_ints[test_bitmap] = wah_to_ints(wah_bms[test_bitmap],
                                                    wah_sizes[test_bitmap],
                                                    &(wah_ints[test_bitmap]));
        }

        wah_int_i = 0;
        wah_bit_i = 0;

        for (ubin_int_i = 0; ubin_int_i < num_ints; ++ubin_int_i) {
            for (ubin_bit_i = 0; ubin_bit_i < 16; ++ubin_bit_i) {
                two_bit = (gt[ubin_int_i] >> (30-(ubin_bit_i*2))) & 3;

                for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap) {
                    bit = (wah_ints[test_bitmap][wah_int_i] >> 
                            (31-wah_bit_i)) & 1;

                    if (test_bitmap == two_bit)
                        TEST_ASSERT_EQUAL(1,bit);
                    else
                        TEST_ASSERT_EQUAL(0,bit);
                }

                wah_bit_i += 1;

                if (wah_bit_i == 32) {
                    wah_int_i += 1;
                    wah_bit_i = 0;
                }


                field_i += 1;
                if (field_i == wf.num_fields)
                    break;
            }
            if (field_i == wf.num_fields)
                break;
        }
        for (test_bitmap = 0; test_bitmap < 4; ++test_bitmap)
            free(wah_bms[test_bitmap]);
    }

    fclose(wf.file);
    free(wf.record_offsets);

    close_hdf5_file(hdf5_f);
    remove(hdf5_file_name);
}
//}}}
#endif

#if 0
//{{{void test_hdf5(void)
//void test_hdf5(void)
{
    uint32_t num_vars = 43;
    uint32_t num_inds = 10;
    char *hdf5_file_name = "../data/10.1e4.var.h5";
    struct hdf5_file hdf5_f = init_hdf5_file(hdf5_file_name,
                                             num_vars,
                                             num_inds);

    uint32_t num_ints = 1 + ((num_inds - 1) / 32);

    uint32_t in[3][num_ints], out[3][num_ints];

    uint32_t i;
    for (i = 0; i < num_ints; ++i){
        in[0][i] = rand();
        in[1][i] = rand();
        in[2][i] = rand();
    }

    write_hdf5_gt(hdf5_f, 0, in[0], "zero");
    write_hdf5_gt(hdf5_f, 1, in[1], "one");
    write_hdf5_gt(hdf5_f, 2, in[2], "two");

    char *md_out;
    read_hdf5_md(hdf5_f, 0, &md_out);
    TEST_ASSERT_EQUAL(0, strcmp("zero", md_out));
    free(md_out);

    read_hdf5_md(hdf5_f, 1, &md_out);
    TEST_ASSERT_EQUAL(0, strcmp("one", md_out));
    free(md_out);

    read_hdf5_md(hdf5_f, 2, &md_out);
    TEST_ASSERT_EQUAL(0, strcmp("two", md_out));
    free(md_out);

    read_hdf5_gt(hdf5_f, 0, out[0]);
    read_hdf5_gt(hdf5_f, 1, out[1]);
    read_hdf5_gt(hdf5_f, 2, out[2]);


    for (i = 0; i < num_ints; ++i) {
        TEST_ASSERT_EQUAL(in[0][i], out[0][i]);
        TEST_ASSERT_EQUAL(in[0][i], out[0][i]);
        TEST_ASSERT_EQUAL(in[0][i], out[0][i]);
    }

    close_hdf5_file(hdf5_f);

    TEST_ASSERT_EQUAL(0, remove(hdf5_file_name));

}
//}}}
#endif

#if 0
//{{{void test_md_bcf_line
//void test_md_bcf_line(void)
{
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);

    int num_samples, num_gts_per_sample;
    char *md;
    uint32_t md_len;

    int r = read_unpack_next_bcf_line(&bcf_f,
                                      &num_samples,
                                      &num_gts_per_sample);
    md_len = md_bcf_line(bcf_f,
                         &md);

    TEST_ASSERT_EQUAL(0, strcmp("1\t0\tV1\tA\tT", md));

    free(md);

    r = read_unpack_next_bcf_line(&bcf_f,
                                  &num_samples,
                                  &num_gts_per_sample);
    md_len = md_bcf_line(bcf_f,
                         &md);

    TEST_ASSERT_EQUAL(0, strcmp("1\t1\tV2\tA\tT", md));

    free(md);
}
//}}}
#endif

#if 0
//{{{void test_append_hdf5_array(void)
//void test_append_hdf5_array(void)
{
    uint32_t num_vars = 43;
    uint32_t num_inds = 10;

    char *hdf5_file_name = "../data/10.1e4.var.h5";
    struct hdf5_file hdf5_f = init_hdf5_file(hdf5_file_name,
                                             num_vars,
                                             num_inds);


    int r = init_r_gt(hdf5_f);

    uint32_t I[5] = {rand(), rand(), rand(), rand(), rand()};
    uint32_t C[5] = {rand()%2, rand()%2, rand()%2, rand()%2, rand()%2};

    set_r_gt(hdf5_f, 1, C[0], I[0]);
    set_r_gt(hdf5_f, 2, C[1], I[1]);
    set_r_gt(hdf5_f, 3, C[2], I[2]);
    set_r_gt(hdf5_f, 4, C[3], I[3]);
    set_r_gt(hdf5_f, 5, C[4], I[4]);

    uint32_t *r_gt_out =
        (uint32_t *)malloc(hdf5_f.num_r_gt_ints * sizeof(uint32_t));

    uint32_t i;

    for (i = 0; i < 5; ++i) {
        read_hdf5_r_gt(hdf5_f, i+1, r_gt_out);
        TEST_ASSERT_EQUAL(I[i], r_gt_out[C[i]]);
    }
    close_hdf5_file(hdf5_f);
}
//}}}
#endif

#if 0
//{{{void test_pack_sum_count_prefix_bcf_line(void)
//void test_pack_sum_count_prefix_bcf_line(void)
{
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);

    int num_samples, num_gts_per_sample;

    int r = read_unpack_next_bcf_line(&bcf_f,
                                      &num_samples,
                                      &num_gts_per_sample);

    uint32_t *packed_ints;
    uint32_t sum, prefix_len, packed_ints_len;

    packed_ints_len = pack_sum_count_prefix_bcf_line(bcf_f,
                                                   num_samples,
                                                   num_gts_per_sample,
                                                   &packed_ints,
                                                   &sum,
                                                   &prefix_len);


    TEST_ASSERT_EQUAL(1, packed_ints_len);
    TEST_ASSERT_EQUAL(2420379648, packed_ints[0]);
    TEST_ASSERT_EQUAL(6, sum);
    TEST_ASSERT_EQUAL(0, prefix_len);

    free(packed_ints);

    r = read_unpack_next_bcf_line(&bcf_f,
                                  &num_samples,
                                  &num_gts_per_sample);

    packed_ints_len = pack_sum_count_prefix_bcf_line(bcf_f,
                                                   num_samples,
                                                   num_gts_per_sample,
                                                   &packed_ints,
                                                   &sum,
                                                   &prefix_len);

    TEST_ASSERT_EQUAL(1, packed_ints_len);
    TEST_ASSERT_EQUAL(16384, packed_ints[0]);
    TEST_ASSERT_EQUAL(1, sum);
    TEST_ASSERT_EQUAL(17, prefix_len);

    free(packed_ints);
}
//}}}
#endif
#if 0
//{{{void test_bcf_read(void)
//void test_bcf_read(void)
{
    char *bcf_file_name = "../data/10.1e4.var.bcf";
    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);

    TEST_ASSERT_EQUAL(10, bcf_f.num_records);

    int num_samples, num_gts_per_sample;

    int r = read_unpack_next_bcf_line(&bcf_f,
                                      &num_samples,
                                      &num_gts_per_sample);

    TEST_ASSERT_EQUAL(10, num_samples);
    TEST_ASSERT_EQUAL(2, num_gts_per_sample);

    int32_t *gt_i = bcf_f.gt;

    TEST_ASSERT_EQUAL(0,
                      strcmp("1",bcf_hdr_id2name(bcf_f.hdr, bcf_f.line->rid)));

    TEST_ASSERT_EQUAL(0, bcf_f.line->pos);
    TEST_ASSERT_EQUAL(0, strcmp("V1", bcf_f.line->d.id));
    TEST_ASSERT_EQUAL(0, strcmp("A", bcf_f.line->d.allele[0]));
    TEST_ASSERT_EQUAL(0, strcmp("T", bcf_f.line->d.allele[1]));

    uint32_t A[20] = { 1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1 };

    int i, j, a = 0;
    for (i = 0; i < num_samples; ++i) {
        for (j=0; j< num_gts_per_sample; ++j) {
            TEST_ASSERT_EQUAL(A[a], bcf_gt_allele(gt_i[j]));
            a+=1;
        }
        gt_i += num_gts_per_sample;
    }
}
//}}}
#endif

//}}}
