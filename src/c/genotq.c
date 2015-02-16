/**
 * @file genotq.c
 * @Author Ryan Layer (ryan.layer@gmail.com)
 * @date May, 2014
 * @brief Functions for converting and opperation on various encoding
 * strategies
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <math.h>
#include <limits.h>
#include "genotq.h"
#include <nmmintrin.h>

// utils
//{{{ void parse_cmd_line_int_csv(int *I,
void parse_cmd_line_int_csv(uint32_t *I,
                            int num_I,
                            char *cmd_line_arg)
{
    char *pch;
    pch = strtok(cmd_line_arg,",");
    int i;
    for (i = 0; i < num_I; ++i){
        I[i] = atoi(pch);
        pch = strtok(NULL,",");
    }
}
//}}}

//{{{uint32_t ints_to_rle(uint32_t *I, int I_len, uint32_t **O)
uint32_t ints_to_rle(uint32_t *I, int I_len, uint32_t **O)
{
    struct uint_ll *head=NULL,*tail=NULL;

    int i,
        j,
        curr_bit,
        ll_len = 0,
        last_bit = -1;
    uint32_t rle_v;
    for (i = 0; i < I_len; ++i) {
        for(j = 0; j < 32; ++j) {
            curr_bit = ((I[i] >> (31-j)) & 1);

            // first one
            if (last_bit == -1) { 
                rle_v = curr_bit << 31;
            // not full, add on
            } else if (curr_bit == last_bit) {
                rle_v += 1;
            // 2^31 -1, this one is full
            } else if ( ((rle_v >> 1) == 2147483647) ||  
                        (curr_bit != last_bit) ) { // diff bit

                struct uint_ll *n = (struct uint_ll *) 
                        malloc(sizeof(struct uint_ll));
                n->value = rle_v;
                n->next = NULL;
                ll_len += 1;

                if (head == NULL)
                    head = n;
                else
                    tail->next = n;

                tail = n;

                rle_v = curr_bit << 31;
            }

            last_bit = curr_bit;
        }
    }

    struct uint_ll *n = (struct uint_ll *) malloc(sizeof(struct uint_ll));
    n->value = rle_v;
    n->next = NULL;
    ll_len += 1;

    if (head == NULL)
        head = n;
    else
        tail->next = n;

    tail = n;

    *O = (uint32_t *) malloc(ll_len*sizeof(uint32_t));
    struct uint_ll *last, *curr = head;
    for (i = 0; i < ll_len; ++i) {
        (*O)[i] = curr->value;
        last = curr;
        curr = curr->next;
        free(last);
    }

    return ll_len;
}
//}}}

//{{{const char *int_to_binary(int x)
const char *int_to_binary(uint32_t x)
{
    static char b[33];
    b[0] = '\0';

    uint32_t z;
    for (z = (uint32_t)(pow(2,31)); z > 0; z >>= 1) {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}
//}}}

//{{{ uint32_t bin_char_to_int(char *bin)
uint32_t bin_char_to_int(char *bin)
{
    uint32_t i = 0;
    int j = 0;

    while (bin[j] != '\0') {
        i = i << 1;
        if (bin[j] == '1')
            i += 1;
        j+=1;
    }

    return i;
}
//}}}

//{{{uint32_t *unpack_1_bit_ints(int packed_int)
int *unpack_1_bit_ints(uint32_t packed_ints)
{
    int *r = (int *) malloc (32*sizeof(int));

    int i;
    for (i = 0; i < 32; ++i) 
        r[i] = (packed_ints >> (31 - i)) & 1;
    
    return r;
}
//}}}

//{{{uint32_t *unpack_2_bit_ints(int packed_int)
int *unpack_2_bit_ints(uint32_t packed_ints)
{
    int *r = (int *) malloc (16*sizeof(int));

    int i;
    for (i = 0; i < 16; ++i) 
        r[i] = (packed_ints >> (30 - i*2)) & 3;
    
    return r;
}
//}}}

//{{{int popcount(uint32_t x) {
int popcount(uint32_t x) {
    /*
    int count;
    for (count=0; x; count++)
        x &= x-1;
    return count;
    */

    return _mm_popcnt_u32(x);

}
//}}}
