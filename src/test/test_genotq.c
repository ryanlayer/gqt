#include <stdio.h>
#include <stdlib.h>
#include "genotq.h"
#include "unity.h"
#include <math.h>

struct plt_file pltf_ind, pltf_var;

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


void setUp(void)
{
    pltf_ind = init_plt_file("data/10.1e4.ind.txt");
    pltf_var = init_plt_file("data/10.1e4.var.txt");
}

void tearDown(void)
{
    fclose(pltf_ind.file);
    fclose(pltf_var.file);
}

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
    TEST_ASSERT_EQUAL(43, pltf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, pltf_ind.num_records);

    TEST_ASSERT_EQUAL(10, pltf_var.num_fields);
    TEST_ASSERT_EQUAL(43, pltf_var.num_records);
}
//}}}

//{{{void test_or_records_vs_field_or(void)
void test_or_records_vs_field_or(void)
{
    int q[4] = {0,1,2,3};
    
    int *R = (int *) calloc(pltf_var.num_fields, sizeof(int));
    int r = or_records_plt(pltf_var, q, 4, R);

    int *F = (int *) calloc(pltf_ind.num_records, sizeof(int));
    int f = or_fields_plt(pltf_ind, q, 4, F);

    int i;
    for (i = 0; i < pltf_ind.num_records; ++i)
        TEST_ASSERT_EQUAL(R[i], F[i]);


    TEST_ASSERT_EQUAL(R[0] , 3);
    TEST_ASSERT_EQUAL(R[1] , 1);
    TEST_ASSERT_EQUAL(R[2] , 0);
    TEST_ASSERT_EQUAL(R[3] , 0);
    TEST_ASSERT_EQUAL(R[4] , 1);
    TEST_ASSERT_EQUAL(R[5] , 0);
    TEST_ASSERT_EQUAL(R[6] , 1);
    TEST_ASSERT_EQUAL(R[7] , 0);
    TEST_ASSERT_EQUAL(R[8] , 1);
    TEST_ASSERT_EQUAL(R[9] , 1);

    free(R);
    free(F);
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
    int len;

    int x = plt_line_to_packed_ints(l1, 43, &r, &len);
    TEST_ASSERT_EQUAL(len, 3);
    TEST_ASSERT_EQUAL(r[0], 2232680464);
    TEST_ASSERT_EQUAL(r[1], 2420396373);
    TEST_ASSERT_EQUAL(r[2], 1163202560);

    free(r);
    
    
    char *l2="1 0 0 0 0 0 1 1 1 0 1 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 1 0 1 1 1 1 1 1 0 1 0 1 1 0 1 0 0";
    //01000000000001010100010101000100 -> 1074087236
    //01000001000100010101010001010101 -> 1091654741
    //01010001000101000100000000000000 -> 1360281600
    
    x = plt_line_to_packed_ints(l2, 43, &r, &len);
    TEST_ASSERT_EQUAL(len, 3);
    TEST_ASSERT_EQUAL(r[0], 1074087236);
    TEST_ASSERT_EQUAL(r[1], 1091654741);
    TEST_ASSERT_EQUAL(r[2], 1360281600);

    free(r);
}
//}}}

//{{{ void test_plot_to_ubin(void) 
void test_plot_to_ubin(void) 
{

    //char *in_file_name="data/10.1e4.ind.txt";
    char *out_file_name="data/.tmp";

    plt_to_ubin(pltf_ind, out_file_name);
}
//}}}

//{{{ void test_init_ubin_file(void)
void test_init_ubin_file(void)
{
    TEST_ASSERT_EQUAL(43, pltf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, pltf_ind.num_records);

    char *out_file_name_1="data/.tmp1";
    plt_to_ubin(pltf_ind, out_file_name_1);

    struct ubin_file uf_ind = init_ubin_file(out_file_name_1);
    TEST_ASSERT_EQUAL(43, uf_ind.num_fields);
    TEST_ASSERT_EQUAL(10, uf_ind.num_records);

    fclose(uf_ind.file);


    TEST_ASSERT_EQUAL(10, pltf_var.num_fields);
    TEST_ASSERT_EQUAL(43, pltf_var.num_records);

    char *out_file_name_2="data/.tmp2";
    plt_to_ubin(pltf_var, out_file_name_2);

    struct ubin_file uf_var = init_ubin_file(out_file_name_2);
    TEST_ASSERT_EQUAL(10, uf_var.num_fields);
    TEST_ASSERT_EQUAL(43, uf_var.num_records);

    fclose(uf_var.file);
}
//}}} 

//{{{ void test_or_records_plt_vs_ubin(void)
void test_or_records_plt_vs_ubin(void)
{
    int q[4] = {0,1,2,3};
    
    int *plt_or = (int *) calloc(pltf_var.num_fields, sizeof(int));
    int plt_or_size = or_records_plt(pltf_var, q, 4, plt_or);

    TEST_ASSERT_EQUAL(plt_or_size, pltf_var.num_fields);

    char *out_file_name_2="data/.tmp2";
    plt_to_ubin(pltf_var, out_file_name_2);

    struct ubin_file uf_var = init_ubin_file(out_file_name_2);
    TEST_ASSERT_EQUAL(10, uf_var.num_fields);
    TEST_ASSERT_EQUAL(43, uf_var.num_records);

    unsigned int *ubin_or;
    int ubin_or_size =  or_records_ubin(uf_var, q, 4, &ubin_or);
    fclose(uf_var.file);

    int *u = unpack_2_bit_ints(ubin_or[0]);

    int i;
    for (i = 0; i < uf_var.num_fields; ++i)
        TEST_ASSERT_EQUAL(u[i], plt_or[i]);

    free(plt_or);
    free(ubin_or);
    free(u);
}
//}}}

//{{{ void test_or_fields_ubin(void)
void test_or_fields_ubin(void)
{
    int q[4] = {0,1,2,3};
    
    int *plt_or = (int *) calloc(pltf_var.num_records, sizeof(int));
    int plt_or_size = or_fields_plt(pltf_var, q, 4, plt_or);

    TEST_ASSERT_EQUAL(plt_or_size, pltf_var.num_records);

    char *out_file_name_2="data/.tmp2";
    plt_to_ubin(pltf_var, out_file_name_2);

    struct ubin_file uf_var = init_ubin_file(out_file_name_2);
    TEST_ASSERT_EQUAL(10, uf_var.num_fields);
    TEST_ASSERT_EQUAL(43, uf_var.num_records);

    unsigned int *ubin_or;
    int ubin_or_size =  or_fields_ubin(uf_var, q, 4, &ubin_or);
    fclose(uf_var.file);

    int *u = unpack_2_bit_ints(ubin_or[0]);

    int i;
    for (i = 0; i < uf_var.num_fields; ++i) {
        TEST_ASSERT_EQUAL(u[i], plt_or[i]);
    }

    free(plt_or);
    free(ubin_or);
    free(u);
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
    unsigned int num_31_groups = map_from_32_bits_to_31_bits(I,5,&O);

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
    
    num_31_groups = map_from_32_bits_to_31_bits(I2,6,&O);

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
    unsigned int wah_size = ints_to_wah(I,5,&O);

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

    wah_size = ints_to_wah(I2,6,&O);

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
    wah_size = ints_to_wah(I3,6,&O);

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
    unsigned int wah_size = ints_to_wah(I,5,&O);

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

    wah_size = ints_to_wah(I2,6,&O);

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
    int wah_size_X = ints_to_wah(X,5,&w_X);
    TEST_ASSERT_EQUAL(5, wah_size_X);
    struct wah_run r_X = init_wah_run(w_X, wah_size_X);

    unsigned int *w_Y;
    int wah_size_Y = ints_to_wah(Y,5,&w_Y);
    TEST_ASSERT_EQUAL(4, wah_size_Y);
    struct wah_run r_Y = init_wah_run(w_Y, wah_size_Y);

    unsigned int *Z;
    unsigned int Z_len = wah_or(&r_X, &r_Y, &Z);

    int i;
    for (i = 0; i < Z_len; ++i)
        TEST_ASSERT_EQUAL(A[i], Z[i]);

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
    unsigned int wah_size = ints_to_wah(I,5,&WAH);

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

    wah_size = ints_to_wah(I2,6,&WAH);
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

    wah_size = ints_to_wah(I3,6,&WAH);
    ints_size = wah_to_ints(WAH,wah_size,&INTS);

    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(I3[i], INTS[i]);

    free(INTS);
    free(WAH);
}
//}}}
