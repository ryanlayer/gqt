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
#include <pthread.h>
#include <assert.h>
#include <immintrin.h>
#include "genotq.h"
#include "pthread_pool.h"

const int tab32[32] = {
    0,  9,  1, 10, 13, 21,  2, 29,
    11, 14, 16, 18, 22, 25,  3, 30,
    8, 12, 20, 28, 15, 17, 24,  7,
    19, 27, 23,  6, 26,  5,  4, 31};

int log2_32 (uint32_t value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return tab32[(uint32_t)(value*0x07C4ACDD) >> 27];
}


// wahbm
//{{{ struct wah_file init_wahbm_file(char *file_name)
struct wah_file init_wahbm_file(char *file_name)
{
    struct wah_file wf;

    wf.file = fopen(file_name, "rb");

    if (!wf.file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        return wf;
    }

    // Jump to the begining of the file to grab the record size
    fseek(wf.file, 0, SEEK_SET);
    fread(&wf.num_fields,sizeof(unsigned int),1,wf.file);
    fread(&wf.num_records,sizeof(unsigned int),1,wf.file);

    wf.record_offsets = (unsigned int *) 
            malloc(sizeof (unsigned int)*wf.num_records*4);

    unsigned int i;
    for (i = 0; i < wf.num_records*4; ++i)
        fread(&(wf.record_offsets[i]),sizeof(unsigned int),1,wf.file);


    wf.header_offset = ftell(wf.file);

    return wf;
}
//}}}

//{{{ unsigned int print_wahbm(struct wah_file wf,
unsigned int print_wahbm(struct wah_file wf,
                         unsigned int *record_ids,
                         unsigned int num_r,
                         unsigned int format)
{
    unsigned int i,j,k,l, bm_size, to_print = num_r;
    unsigned int *bm = NULL;

    unsigned int num_ints_per_record = 1 + ((wf.num_fields - 1) / 16);

    if (num_r == 0)
        to_print = wf.num_records;


    unsigned int *output_record = (unsigned int *) malloc 
            (num_ints_per_record * sizeof(unsigned int));

    unsigned int *tmp_record = (unsigned int *) malloc 
            (num_ints_per_record * sizeof(unsigned int));

    for (i = 0; i < to_print; ++i) {

        memset(output_record, 0, num_ints_per_record * sizeof(unsigned int));

        for (j = 0; j < 4; ++j) {
            memset(tmp_record, 0, 
                    num_ints_per_record * sizeof(unsigned int));

            // get the compressed bitmap
            if (num_r > 0)
                bm_size = get_wah_bitmap(wf,
                                         record_ids[i],
                                         j,
                                         &bm);
            else
                bm_size = get_wah_bitmap(wf,
                                         i,
                                         j,
                                         &bm);

            // decompress 
            unsigned int *ints = NULL;
            unsigned int ints_size = wah_to_ints(bm,bm_size,&ints);


#if 0

            for (k = 0; k < ints_size; ++k) {
                if (k !=0)
                    printf(" ");
                for (l = 0; l < 32; ++l) {
                    if (l !=0)
                        printf(" ");
                    unsigned int bit = (ints[k] >> (31 - l)) & 1;
                    printf("%u",bit);
                    
                }
            }
            printf("\n");
#endif



#if 1
            // loop through each bit, and set the corresponding possition to j
            // if the bit is one
            int int_i = 0, bit_i = 0;
            for (k = 0; k < ints_size; ++k) {
                for (l = 0; l < 32; ++l) {
                    unsigned int bit = (ints[k] >> (31 - l)) & 1;

                    if (bit == 1)
                        tmp_record[int_i] += j << (30 - (bit_i * 2));

                    bit_i += 1;
                    if (bit_i == 16) {
                        int_i += 1;
                        bit_i = 0;
                    }

                }
            }
#endif

            free(bm);
            free(ints);
            bm = NULL;
            ints = NULL;

#if 1
            for (k = 0; k < num_ints_per_record; ++k) 
                output_record[k] += tmp_record[k];
#endif
        }

#if 1
        unsigned int printed_bits = 0;
        for (j = 0; j < num_ints_per_record; ++j) {
            if (j !=0)
                printf(" ");
            for (k = 0; k < 16; ++k) {
                unsigned int bit = (output_record[j] >> (30 - 2*k)) & 3;
                if (k !=0)
                    printf(" ");
                printf("%u", bit);
                printed_bits += 1;
                if (printed_bits == wf.num_fields)
                    break;
            }
        }
        printf("\n");
#endif
    }

    free(tmp_record);
    free(output_record);

    return to_print;
}
//}}}

//{{{ unsigned int print_by_name_wahbm(char *wahbm_file_name,
unsigned int print_by_name_wahbm(char *wahbm_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format)
{
    struct wah_file wf = init_wahbm_file(wahbm_file_name);
    return print_wahbm(wf, record_ids, num_r, format);
}
//}}} 

//{{{ unsigned int get_wah_bitmap(struct wah_file wf,
unsigned int get_wah_bitmap(struct wah_file wf,
                            unsigned int wah_record,
                            unsigned int bitmap,
                            unsigned int **wah_bitmap)
{
    // get the size of the WAH-encoded bitmap
    unsigned int wah_size = 0, wah_offset = 0;
    if ((wah_record == 0) && (bitmap == 0)) {
        wah_size = wf.record_offsets[wah_record + bitmap];
        wah_offset = wf.header_offset;
    } else {
        wah_size = wf.record_offsets[wah_record*4 + bitmap] - 
                   wf.record_offsets[wah_record*4 + bitmap - 1];
        /*
        fprintf(stderr, "wf.header_offset:%lu\t"
                        "wah_record:%u\t"
                        "bitmap:%u\t"
                        "wah_size:%u\t"
                        "wf.record_offsets[]:%u\n",
                        wf.header_offset,
                        wah_record,
                        bitmap,
                        wah_size,
                        wf.record_offsets[wah_record*4 + bitmap]);
        */

        wah_offset = wf.header_offset +
                     sizeof(unsigned int) * 
                        (wf.record_offsets[wah_record*4 + bitmap] - wah_size);
    }

    //fprintf(stderr, "wah_size:%u\twah_offset:%u\n", wah_size, wah_offset);


    *wah_bitmap = (unsigned int *) malloc(sizeof(unsigned int)*wah_size);
    fseek(wf.file, wah_offset, SEEK_SET);
    fread(*wah_bitmap,sizeof(unsigned int),wah_size,wf.file);

    return wah_size;
}
//}}}

//{{{ unsigned int range_records_w_exclude_wahbm(struct wah_file wf,
unsigned int range_records_w_exclude_wahbm(struct wah_file wf,
                                           unsigned int *record_ids,
                                           unsigned int num_r,
                                           unsigned int start_test_value,
                                           unsigned int end_test_value,
                                           unsigned int exclude_value,
                                           unsigned int **R) 

{
    unsigned int *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    unsigned int record_curr_bm_size,
                 record_new_bm_size,
                 record_tmp_bm_size;

    unsigned int *query_curr_bm = NULL,
                 *query_tmp_bm = NULL;

    unsigned int query_curr_bm_size,
                 query_tmp_bm_size;


    unsigned int i,j,k,l;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = start_test_value; j < end_test_value; ++j) {

            if (j == exclude_value)
            {
                continue;
            }

            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);

            if (record_curr_bm == NULL) {
                record_curr_bm = record_new_bm;
                record_curr_bm_size = record_new_bm_size;
            } else {
                struct wah_run curr_run = init_wah_run(record_curr_bm,
                                                       record_curr_bm_size);
                struct wah_run new_run = init_wah_run(record_new_bm,
                                                      record_new_bm_size);

                record_tmp_bm_size = wah_or(&curr_run,
                                            &new_run,
                                            &record_tmp_bm);
                free(record_curr_bm);
                free(record_new_bm);

                record_curr_bm = record_tmp_bm;
                record_curr_bm_size = record_tmp_bm_size;
            }
        }

        if (query_curr_bm == NULL) {
            query_curr_bm = record_curr_bm;
            query_curr_bm_size = record_curr_bm_size;
        } else {
                struct wah_run record_run = init_wah_run(record_curr_bm,
                                                         record_curr_bm_size);
                struct wah_run query_run = init_wah_run(query_curr_bm,
                                                        query_curr_bm_size);

                query_tmp_bm_size = wah_and(&record_run,
                                            &query_run,
                                            &query_tmp_bm);
                free(record_curr_bm);
                free(query_curr_bm);

                query_curr_bm = query_tmp_bm;
                query_curr_bm_size = query_tmp_bm_size;
        }
    }

    *R = query_curr_bm;
    return query_curr_bm_size;
}
//}}}

//{{{ unsigned int count_range_records_wahbm(struct wah_file wf,
unsigned int count_range_records_wahbm(struct wah_file wf,
                                       unsigned int *record_ids,
                                       unsigned int num_r,
                                       unsigned int start_test_value,
                                       unsigned int end_test_value,
                                       unsigned int **R) 

{
    *R = (unsigned int *) calloc(wf.num_fields,sizeof(unsigned int));

    unsigned int *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    unsigned int record_curr_bm_size,
                 record_new_bm_size,
                 record_tmp_bm_size;

    unsigned int i,j, r_size;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = start_test_value; j < end_test_value; ++j) {

            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);

            if (record_curr_bm == NULL) {
                record_curr_bm = record_new_bm;
                record_curr_bm_size = record_new_bm_size;
            } else {
                struct wah_run curr_run = init_wah_run(record_curr_bm,
                                                       record_curr_bm_size);
                struct wah_run new_run = init_wah_run(record_new_bm,
                                                      record_new_bm_size);

                record_tmp_bm_size = wah_or(&curr_run,
                                            &new_run,
                                            &record_tmp_bm);
                free(record_curr_bm);
                free(record_new_bm);

                record_curr_bm = record_tmp_bm;
                record_curr_bm_size = record_tmp_bm_size;
            }
        }

        r_size = add_wahbm(*R,
                           wf.num_fields,
                           record_curr_bm,
                           record_curr_bm_size);
    }
    return wf.num_fields;
}
//}}}

//{{{ unsigned int sum_range_records_wahbm(struct wah_file wf,
unsigned int sum_range_records_wahbm(struct wah_file wf,
                                     unsigned int *record_ids,
                                     unsigned int num_r,
                                     unsigned int **R) 

{
    *R = (unsigned int *) calloc(wf.num_fields,sizeof(unsigned int));

    unsigned int *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    unsigned int record_curr_bm_size,
                 record_new_bm_size,
                 record_tmp_bm_size;

    unsigned int i,j, r_size;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = 0; j < 4; ++j) {

            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);

            if (record_curr_bm == NULL) {
                record_curr_bm = record_new_bm;
                record_curr_bm_size = record_new_bm_size;
            } else {
                struct wah_run curr_run = init_wah_run(record_curr_bm,
                                                       record_curr_bm_size);
                struct wah_run new_run = init_wah_run(record_new_bm,
                                                      record_new_bm_size);

                record_tmp_bm_size = wah_or(&curr_run,
                                            &new_run,
                                            &record_tmp_bm);
                free(record_curr_bm);
                free(record_new_bm);

                record_curr_bm = record_tmp_bm;
                record_curr_bm_size = record_tmp_bm_size;
            }
        }

        r_size = add_wahbm(*R,
                           wf.num_fields,
                           record_curr_bm,
                           record_curr_bm_size);
    }
    return wf.num_fields;
}
//}}}

//{{{ unsigned int range_records_wahbm(struct wah_file wf,
unsigned int range_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int **R) 

{
    unsigned int *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    unsigned int record_curr_bm_size,
                 record_new_bm_size,
                 record_tmp_bm_size;

    unsigned int *query_curr_bm = NULL,
                 *query_tmp_bm = NULL;

    unsigned int query_curr_bm_size,
                 query_tmp_bm_size;


    unsigned int i,j,k,l;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = start_test_value; j < end_test_value; ++j) {

            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);

            if (record_curr_bm == NULL) {
                record_curr_bm = record_new_bm;
                record_curr_bm_size = record_new_bm_size;
            } else {
                struct wah_run curr_run = init_wah_run(record_curr_bm,
                                                       record_curr_bm_size);
                struct wah_run new_run = init_wah_run(record_new_bm,
                                                      record_new_bm_size);

                record_tmp_bm_size = wah_or(&curr_run,
                                            &new_run,
                                            &record_tmp_bm);
                free(record_curr_bm);
                free(record_new_bm);

                record_curr_bm = record_tmp_bm;
                record_curr_bm_size = record_tmp_bm_size;
            }
        }

        if (query_curr_bm == NULL) {
            query_curr_bm = record_curr_bm;
            query_curr_bm_size = record_curr_bm_size;
        } else {
                struct wah_run record_run = init_wah_run(record_curr_bm,
                                                         record_curr_bm_size);
                struct wah_run query_run = init_wah_run(query_curr_bm,
                                                        query_curr_bm_size);

                query_tmp_bm_size = wah_and(&record_run,
                                            &query_run,
                                            &query_tmp_bm);
                free(record_curr_bm);
                free(query_curr_bm);

                query_curr_bm = query_tmp_bm;
                query_curr_bm_size = query_tmp_bm_size;
        }
    }

    *R = query_curr_bm;
    return query_curr_bm_size;
}
//}}}

//{{{ unsigned int add_wahbm(unsigned int *R,
unsigned int add_wahbm(unsigned int *R,
                       unsigned int r_size,
                       unsigned int *wah,
                       unsigned int wah_size)
{

    unsigned int wah_c,
                 wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    field_i = 0;

    unsigned int v;

    for (wah_i = 0; wah_i < wah_size; ++wah_i) {
        wah_c = wah[wah_i];
        if (wah_c >> 31 == 1) {
            num_words = (wah_c & 0x3fffffff);
            fill_bit = (wah_c>=0xC0000000?1:0);
            bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = wah_c;
        }

        if ( (num_words > 1) && (fill_bit == 0) ) {
            field_i += num_words;
            if (field_i >= r_size)
                return r_size;
        } else {
            if (bits == 0) {
                field_i += 31;
                if (field_i >= r_size)
                    return r_size;
            } else {
                for (word_i = 0; word_i < num_words; ++word_i) {
                    /* 
                    // Attempt to reduce the number of times the for loop
                    // itterates so that 
                    v = bits;
                    for ( ; v ; ) {
                        R[field_i] += log2_32(v&(-v));
                        v &= v - 1;
                        field_i += 1;
                        if (field_i >= r_size)
                            return r_size;
                    }
                    */
                    for (bit_i = 0; bit_i < 31; ++bit_i) {
                        R[field_i] += (bits >> (30 - bit_i)) & 1;
                        field_i += 1;

                        if (field_i >= r_size)
                            return r_size;
                    }
                }
            }
        }
    }

    return r_size;
}
//}}}

#ifdef __AVX2__
//{{{void avx_add(unsigned int bits,
void avx_add(unsigned int bits,
             __m256i *s_1,
             __m256i *s_2,
             __m256i *s_3,
             __m256i *s_4,
             __m256i *m,
             __m256i *R_avx,
             unsigned int field_i)
{

    unsigned int avx_i = field_i/8;

    __m256i y1 = _mm256_set1_epi32(bits);

    __m256i y2 = _mm256_srlv_epi32 (y1, *s_1);
    __m256i y3 = _mm256_and_si256 (y2, *m);
    R_avx[3+avx_i] = _mm256_add_epi32(R_avx[3+avx_i], y3);

    y2 = _mm256_srlv_epi32 (y1, *s_2);
    y3 = _mm256_and_si256 (y2, *m);
    R_avx[2+avx_i] = _mm256_add_epi32(R_avx[2+avx_i], y3);

    y2 = _mm256_srlv_epi32 (y1, *s_3);
    y3 = _mm256_and_si256 (y2, *m);
    R_avx[1+avx_i] = _mm256_add_epi32(R_avx[1+avx_i], y3);

    y2 = _mm256_srlv_epi32 (y1, *s_4);
    y3 = _mm256_and_si256 (y2, *m);
    R_avx[0+avx_i] = _mm256_add_epi32(R_avx[0+avx_i], y3);
}
//}}}
#endif

#ifdef __AVX2__
//{{{ unsigned int avx_add_wahbm(unsigned int *R,
unsigned int avx_add_wahbm(unsigned int *R,
                       unsigned int r_size,
                       unsigned int *wah,
                       unsigned int wah_size)
{
    __attribute__((aligned(64))) int rshift_4[8] = { 31, 30, 29, 28, 27, 26, 25, 24 };
    __attribute__((aligned(64))) int rshift_3[8] = { 23, 22, 21, 20, 19, 18, 17, 16 };
    __attribute__((aligned(64))) int rshift_2[8] = { 15, 14, 13, 12, 11, 10, 9, 8 };
    __attribute__((aligned(64))) int rshift_1[8] = { 7, 6, 5, 4, 3, 2, 1, 0 };
    __attribute__((aligned(64))) int masks[8] =  { 1, 1, 1, 1, 1, 1, 1, 1 };

    
    __m256i *R_avx = (__m256i *)R;

    __m256i *rshift_1_avx = (__m256i *)rshift_1;
    __m256i *rshift_2_avx = (__m256i *)rshift_2;
    __m256i *rshift_3_avx = (__m256i *)rshift_3;
    __m256i *rshift_4_avx = (__m256i *)rshift_4;
    __m256i *masks_avx = (__m256i *)masks;

    __m256i s_1 = _mm256_load_si256(rshift_1_avx);
    __m256i s_2 = _mm256_load_si256(rshift_2_avx);
    __m256i s_3 = _mm256_load_si256(rshift_3_avx);
    __m256i s_4 = _mm256_load_si256(rshift_4_avx);
    __m256i m = _mm256_load_si256(masks_avx);
    __m256i y1, y2, y3;

    unsigned int wah_c,
                 wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    field_i = 0;

    unsigned int buf, buf_empty_bits = 32;

    for (wah_i = 0; wah_i < wah_size; ++wah_i) {
        wah_c = wah[wah_i];
        if (wah_c >> 31 == 1) {
            num_words = (wah_c & 0x3fffffff);
            fill_bit = (wah_c>=0xC0000000?1:0);
            bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = wah_c;
        }

        if ( (num_words > 1) && (fill_bit == 0) ) {
            // probably need to account for extra bits here
            if (buf_empty_bits < 32)
                avx_add(buf, &s_1, &s_2, &s_3, &s_4, &m, R_avx, field_i);
            
            field_i += 32;
            // the empty bits were supplied by this run, so we don't want to
            // count them twice
            field_i += num_words*31 - buf_empty_bits; 

            buf_empty_bits = 32;
            buf = 0;

            if (field_i >= r_size)
                return r_size;
        } else {
            if (bits == 0) {

                if (buf_empty_bits < 32)
                    avx_add(buf, &s_1, &s_2, &s_3, &s_4, &m, R_avx, field_i);
                field_i += 32 + (31 - buf_empty_bits);

                buf = 0;
                buf_empty_bits = 32;

                if (field_i >= r_size)
                    return r_size;
/*
                
                if (buf_empty_bits < 32)
                    avx_add(buf, &s_1, &s_2, &s_3, &s_4, &m, R_avx, field_i);
                field_i += 32;

                buf = 0;
                buf_empty_bits = 32 - buf_empty_bits;
*/

                if (field_i >= r_size)
                    return r_size;

            } else {
                for (word_i = 0; word_i < num_words; ++word_i) {
                    if (buf_empty_bits == 32) {
                        if (field_i % 32 != 0) {
                            // add padding to buf that makes up for the
                            // difference, then add 32 - (field_i % 32) bits to
                            // the buff
                            unsigned int padding = field_i % 32;
                            buf = bits >> (padding - 1);
                            avx_add(buf,
                                    &s_1,
                                    &s_2,
                                    &s_3,
                                    &s_4,
                                    &m,
                                    R_avx,
                                    field_i - padding);
                            field_i+= 32 - padding;
                            buf = bits << (32 - padding + 1);
                            buf_empty_bits = (32 - padding) + 1;
                        } else {
                            buf = bits << 1;
                            buf_empty_bits = 1;
                        }
                    } else {

                        buf += bits >> (31-buf_empty_bits);

                        avx_add(buf,
                                &s_1,
                                &s_2,
                                &s_3,
                                &s_4,
                                &m,
                                R_avx,
                                field_i);

                        field_i+=32;

                        buf_empty_bits += 1;
                        buf = bits << buf_empty_bits;
                    }
                }
            }
        }
    }

    for (bit_i = 0; bit_i < 31; ++bit_i) {
        R[field_i] += (buf >> (31 - bit_i)) & 1;
        field_i += 1;

        if (field_i >= r_size)
            return r_size;
    }

    return r_size;
}
//}}}
#endif

//{{{ unsigned int add_n_wahbm(unsigned int *R,
unsigned int add_n_wahbm(unsigned int *R,
                       unsigned int n,
                       unsigned int r_size,
                       unsigned int *wah,
                       unsigned int wah_size)
{

    unsigned int wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    unsigned int t;
    field_i = 0;

    for (wah_i = 0; wah_i < wah_size; ++wah_i) {

        if (wah[wah_i] >> 31 == 1) {
            num_words = (wah[wah_i] & 0x3fffffff);
            fill_bit = (wah[wah_i]>=0xC0000000?1:0);
            bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = wah[wah_i];
        }

        if ( (num_words > 1) && (fill_bit == 0) ) {
            field_i += num_words;
            if (field_i >= r_size)
                return r_size;
        } else {
            if (bits == 0) {
                field_i += 31;
                if (field_i >= r_size)
                    return r_size;
            } else {

                for (word_i = 0; word_i < num_words; ++word_i) {
#if 0
                    for( ; bits; ) {
                        t = log2_32(bits&(-bits));
                        R[field_i + 30 - t]+=1;
                        bits &= bits - 1;
                    }
                    field_i += 31;
                    if (field_i >= r_size)
                        return r_size;
#endif
                        
#if 1
                    for (bit_i = 0; bit_i < 31; ++bit_i) {
                        bit = (bits >> (30 - bit_i)) & 1;
                        R[field_i] += bit * n;
                        field_i += 1;

                        if (field_i >= r_size)
                            return r_size;
                    }
#endif
                }
            }
        }
    }

    return r_size;
}
//}}}

//{{{void *t_add_n_wahbm(void *arg)
void *t_add_n_wahbm(void *arg)
{
    struct t_add_n_wahbm_args *a = (struct t_add_n_wahbm_args *)arg;

    unsigned int bit_i,
                 bit,
                 bits = a->bits,
                 field_i = a->field_i, 
                 n = a->n, 
                 r_size = a->r_size;

    //fprintf(stderr, "field_i:%u\n", field_i);
    for (bit_i = 0; bit_i < 31; ++bit_i) {
        bit = (bits >> (30 - bit_i)) & 1;
        a->R[field_i] += bit * n;
        field_i += 1;

        if (field_i >= r_size)
            return NULL;
    }

    return NULL;
}
//}}}

//{{{void *t_add_n_wahbm_2(void *arg)
void *t_add_n_wahbm_2(void *arg)
{
    struct t_add_n_wahbm_args *a = (struct t_add_n_wahbm_args *)arg;

    unsigned int bit_i,
                 word_i,
                 bit,
                 bits = a->bits,
                 field_i = a->field_i, 
                 n = a->n, 
                 r_size = a->r_size,
                 num_words = a->num_words;

    for (word_i = 0; word_i < num_words; ++word_i) {
        //fprintf(stderr, "field_i:%u\n", field_i);
        for (bit_i = 0; bit_i < 31; ++bit_i) {
            bit = (bits >> (30 - bit_i)) & 1;
            a->R[field_i] += bit * n;
            field_i += 1;

            if (field_i >= r_size)
                return NULL;
        }
    }

    return NULL;
}
//}}}

//{{{ unsigned int p_pool_add_n_wahbm(unsigned int *R,
unsigned int p_pool_add_n_wahbm(unsigned int *R,
                                unsigned int n,
                                unsigned int r_size,
                                unsigned int *wah,
                                unsigned int wah_size,
                                struct pool *t_pool)
{
    unsigned int wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    unsigned int t;
    field_i = 0;

    struct t_add_n_wahbm_args *arg;

    for (wah_i = 0; wah_i < wah_size; ++wah_i) {

        if (wah[wah_i] >> 31 == 1) {
            num_words = (wah[wah_i] & 0x3fffffff);
            fill_bit = (wah[wah_i]>=0xC0000000?1:0);
            bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = wah[wah_i];
        }

        if ( (num_words > 1) && (fill_bit == 0) ) {
            field_i += num_words;
            if (field_i >= r_size)
                return r_size;
        } else {
            if (bits == 0) {
                field_i += 31;
                if (field_i >= r_size)
                    return r_size;
            } else {

#if 1

                    arg = (struct t_add_n_wahbm_args *)
                            malloc(sizeof(struct t_add_n_wahbm_args));
                    arg->bits = bits;
                    arg->field_i = field_i;
                    arg->R = R;
                    arg->r_size = r_size;
                    arg->n = n;
                    arg->n = num_words;

                    pool_enqueue(t_pool, arg, 1);

                    field_i += num_words*31;
                    if (field_i >= r_size)
                        return r_size;
#endif

#if 0
                for (word_i = 0; word_i < num_words; ++word_i) {

                    arg = (struct t_add_n_wahbm_args *)
                            malloc(sizeof(struct t_add_n_wahbm_args));
                    arg->bits = bits;
                    arg->field_i = field_i;
                    arg->R = R;
                    arg->r_size = r_size;
                    arg->n = n;

                    //fprintf(stderr, "Q:");
                    pool_enqueue(t_pool, arg, 1);
                    //fprintf(stderr, "%u\n", field_i);
                    /*
                    for (bit_i = 0; bit_i < 31; ++bit_i) {
                        bit = (bits >> (30 - bit_i)) & 1;
                        R[field_i] += bit * n;
                        field_i += 1;

                        if (field_i >= r_size)
                            return r_size;
                    }
                    */

                    field_i += 31;
                    if (field_i >= r_size)
                        return r_size;
                }
#endif
            }
        }
    }

    return r_size;
}
//}}}

//{{{ eq, ne, gt, gte, lt, lte :records_wahbm
//{{{ unsigned int eq_records_wahbm(struct wah_file wf,
unsigned int eq_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value, test_value+1, R);
}
//}}}

//{{{ unsigned int ne_records_wahbm(struct wah_file wf,
unsigned int ne_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for lower bound and upper bound.
    // exclude the test_value
    return range_records_w_exclude_wahbm(wf, record_ids, num_r, 0, 4, test_value, R);
}
//}}}

//{{{ unsigned int gt_records_wahbm(struct wah_file wf,
unsigned int gt_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value+1, 4, R);
}
//}}}

//{{{ unsigned int gte_records_wahbm(struct wah_file wf,
unsigned int gte_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value, 4, R);
}
//}}}

//{{{ unsigned int lt_records_wahbm(struct wah_file wf,
unsigned int lt_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, 0, test_value, R);
}
//}}}

//{{{ unsigned int lte_records_wahbm(struct wah_file wf,
unsigned int lte_records_wahbm(struct wah_file wf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int test_value,
                              unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, 0, test_value+1, R);
}
//}}}
//}}}

//{{{ unsigned int gt_count_records_wahbm(struct wah_file wf,
unsigned int gt_count_records_wahbm(struct wah_file wf,
                                    unsigned int *record_ids,
                                    unsigned int num_r,
                                    unsigned int test_value,
                                    unsigned int **R)
{
    // TODO: need constants for upper bound.
    return count_range_records_wahbm(wf,
                                     record_ids,
                                     num_r,
                                     test_value+1,
                                     4,
                                     R);
}
//}}}
