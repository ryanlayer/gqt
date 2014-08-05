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


//{{{ struct ubin_file init_ubin_file(char *file_name)
struct ubin_file init_ubin_file(char *file_name)
{
    struct ubin_file uf;

    uf.file = fopen(file_name, "rb");

    if (!uf.file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        return uf;
    }

    // Jump to the begining of the file to grab the record size
    fseek(uf.file, 0, SEEK_SET);
    fread(&uf.num_fields,sizeof(unsigned int),1,uf.file);
    fread(&uf.num_records,sizeof(unsigned int),1,uf.file);
    uf.header_offset = ftell(uf.file);

    return uf;

}
//}}}

//{{{ unsigned int convert_file_by_name_ubin_to_wahbm16(char *ubin_in, 
unsigned int convert_file_by_name_ubin_to_wahbm16(char *ubin_in,
                                                  char *wah_out)
{
    FILE *wf = fopen(wah_out,"wb");

    if (!wf) {
        printf("Unable to open %s\n", wah_out);
        return 1;
    }

    struct ubin_file uf = init_ubin_file(ubin_in);

    //write header for WAH bitmap index file
    fwrite(&(uf.num_fields), sizeof(int), 1, wf);
    fwrite(&(uf.num_records), sizeof(int), 1, wf);
    int zero = 0;
    int k;
    for (k = 0; k < uf.num_records*4; ++k)
        fwrite(&zero, sizeof(int), 1, wf);

    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    int i,j,wah_i = 0, offset_total  = 0;

    // skip to the target record and read in the full record
    fseek(uf.file, uf.header_offset, SEEK_SET);

    for (i = 0; i < uf.num_records; ++i) {
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);
         
        uint16_t *wah;
        unsigned int *wah_sizes;
        unsigned int wah_len = ubin_to_bitmap_wah16(c,
                                                    num_ints_per_record,
                                                    uf.num_fields,
                                                    &wah,
                                                    &wah_sizes);

        fseek(wf,sizeof(unsigned int)* (2+4* wah_i),  SEEK_SET);
        for (j = 0; j < 4; ++j) {
            offset_total += wah_sizes[j];
            fwrite(&offset_total, sizeof(unsigned int), 1, wf);
        }

        fseek(wf,0,SEEK_END);
        size_t ret = fwrite(wah, sizeof(uint16_t), wah_len, wf);
        if (ret != wah_len)
            fprintf(stderr, "ret:%zu != wah_len:%u\n", ret, wah_len);

        wah_i+=1;
        free(wah);
        free(wah_sizes);
    }

    free(c);

    fclose(wf);
    fclose(uf.file);
    return 0;
}
//}}}

//{{{ unsigned int convert_file_by_name_ubin_to_wahbm(char *ubin_in, 
unsigned int convert_file_by_name_ubin_to_wahbm(char *ubin_in, char *wah_out)
{
    FILE *wf = fopen(wah_out,"wb");

    if (!wf) {
        printf("Unable to open %s\n", wah_out);
        return 1;
    }

    struct ubin_file uf = init_ubin_file(ubin_in);

    //write header for WAH bitmap index file
    fwrite(&(uf.num_fields), sizeof(int), 1, wf);
    fwrite(&(uf.num_records), sizeof(int), 1, wf);
    int zero = 0;
    int k;
    for (k = 0; k < uf.num_records*4; ++k)
        fwrite(&zero, sizeof(int), 1, wf);

    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    int i,j,wah_i = 0, offset_total  = 0;

    // skip to the target record and read in the full record
    fseek(uf.file, uf.header_offset, SEEK_SET);

    for (i = 0; i < uf.num_records; ++i) {
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);
         
        unsigned int *wah;
        unsigned int *wah_sizes;
        unsigned int wah_len = ubin_to_bitmap_wah(c,
                                                  num_ints_per_record,
                                                  uf.num_fields,
                                                  &wah,
                                                  &wah_sizes);

        fseek(wf,sizeof(unsigned int)* (2+4* wah_i),  SEEK_SET);
        for (j = 0; j < 4; ++j) {
            offset_total += wah_sizes[j];
            fwrite(&offset_total, sizeof(unsigned int), 1, wf);
            //fprintf(stderr,"%u\t%u\n", wah_sizes[j], offset_total);
        }

        fseek(wf,0,SEEK_END);
        size_t ret = fwrite(wah, sizeof(unsigned int), wah_len, wf);
        if (ret != wah_len)
            fprintf(stderr, "ret:%zu != wah_len:%u\n", ret, wah_len);

        wah_i+=1;
        free(wah);
        free(wah_sizes);
    }

    free(c);

    fclose(wf);
    fclose(uf.file);
    return 0;
}
//}}}

//{{{ unsigned int convert_file_by_name_ubin_to_wah(char *ubin_in,
unsigned int convert_file_by_name_ubin_to_wah(char *ubin_in, char *wah_out)
{
    FILE *wf = fopen(wah_out,"wb");

    if (!wf) {
        printf("Unable to open %s\n", wah_out);
        return 1;
    }

    struct ubin_file uf = init_ubin_file(ubin_in);

    //write header for WAH bitmap index file
    fwrite(&(uf.num_fields), sizeof(int), 1, wf);
    fwrite(&(uf.num_records), sizeof(int), 1, wf);
    int zero = 0;
    int k;
    for (k = 0; k < uf.num_records; ++k)
        fwrite(&zero, sizeof(int), 1, wf);

    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    int i,j,wah_i = 0, offset_total  = 0;

    // skip to the target record and read in the full record
    fseek(uf.file, uf.header_offset, SEEK_SET);

    for (i = 0; i < uf.num_records; ++i) {
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);
         
        unsigned int *wah;
        unsigned int wah_len = ints_to_wah(c,
                                           num_ints_per_record,
                                           uf.num_fields*2,
                                           &wah);

        fseek(wf,sizeof(unsigned int)* (2+wah_i),  SEEK_SET);

        offset_total += wah_len;
        fwrite(&offset_total, sizeof(unsigned int), 1, wf);

        fseek(wf,0,SEEK_END);
        size_t ret = fwrite(wah, sizeof(unsigned int), wah_len, wf);
        if (ret != wah_len)
            fprintf(stderr, "ret:%zu != wah_len:%u\n", ret, wah_len);

        wah_i+=1;
        free(wah);
    }

    free(c);

    fclose(wf);
    fclose(uf.file);
    return 0;
}
//}}}

//{{{ unsigned int get_ubin_record(struct ubin_file uf,
unsigned int get_ubin_record(struct ubin_file uf,
                             unsigned int record_id,
                             unsigned int **ubin_record)
{
    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);

    unsigned int ubin_offset = uf.header_offset + 
            sizeof(unsigned int)*(record_id*num_ints_per_record);

    //fprintf(stderr, "ubin_offset:%u\n", ubin_offset);

    *ubin_record = (unsigned int *)
                   malloc(sizeof(unsigned int)*num_ints_per_record);

    fseek(uf.file, ubin_offset, SEEK_SET);
    fread(*ubin_record,sizeof(unsigned int),num_ints_per_record,uf.file);

    return num_ints_per_record;
}

//}}}

//{{{ unsigned int ubin_to_bitmap(unsigned int *U,
unsigned int ubin_to_bitmap(unsigned int *U,
                            unsigned int U_len,
                            unsigned int used_bits,
                            unsigned int **B)
{
    // Since U encodeds a series of two-bit values, and the bitmap uses one
    // bit per unique value in U, the bitmap for each value will require 1/2 
    // (rounded up) the number of ints used to encode U
    unsigned int value_index_size = (U_len + 2 - 1) / 2;
    // There are 4 unique values, so in total B will require 4x  
    unsigned int B_len = 4 * value_index_size;
    *B = (unsigned int *) calloc(B_len, sizeof(unsigned int));

    unsigned int two_bit, set_bit, bit_offset, B_int_i, B_two_bit_i, B_bit_i;

    B_int_i = 0;
    B_two_bit_i = 0;
    B_bit_i = 0;


    unsigned int i,j,k;
    for (i = 0; i < U_len; ++i) {
        for (j = 0; j < 16; ++j) {
            two_bit = ((U[i] >> (30 - (2*j)))& 3);
            for (k = 0; k < 4; ++k) {

                if (k == two_bit)
                    set_bit = 1;
                else 
                    set_bit = 0;

                /* 
                 * B consists of 4 bit arrays laid out consecutively.
                 * As we loop over the two bit values, a single bit must be
                 * set in all 4 bit arrays.  The first bit array is at
                 * B_int_i, and each bit array has value_index_size values
                 *
                 */

                bit_offset = (k * value_index_size) + B_int_i;

                (*B)[bit_offset] +=  set_bit << (31 - B_two_bit_i);
                

                /*
                fprintf(stderr, "k:%u\t"
                                "set_bit:%u\t"
                                "two_bit:%u\t" 
                                "B_int_i:%u\t" 
                                "B_two_bit_i:%u\n",
                                k,
                                set_bit,
                                two_bit,
                                B_int_i,
                                B_two_bit_i);
                */
            }
            B_two_bit_i += 1;
            B_bit_i += 2;

            if (B_two_bit_i == 32) {
                B_int_i += 1;
                B_two_bit_i = 0;
            }

            if (B_bit_i >= used_bits)
                break;
        }
        if (B_bit_i >= used_bits)
            break;
    }

    return B_len;
}
//}}}

//{{{ unsigned int ubin_to_bitmap_wah16(unsigned int *U,
unsigned int ubin_to_bitmap_wah16(unsigned int *U,
                                  unsigned int U_len,
                                  unsigned int num_fields,
                                  uint16_t **W,
                                  unsigned int **wah_sizes)
{
    unsigned int *B;
    // two bits per field
    unsigned int B_len = ubin_to_bitmap(U, U_len, num_fields*2, &B);
    unsigned int b_len = B_len / 4;  // size of each bitmap index

    *wah_sizes = (unsigned int *) malloc(4*sizeof(unsigned int));

    uint16_t *wahs[4];
    unsigned int wahs_size[4],
                 i,
                 j;

    unsigned int total_wah_size = 0;

    for (i = 0; i < 4; i++) {
        wahs_size[i] = ints_to_wah16( (B + (i*b_len)), 
                                      b_len, 
                                      num_fields, // 1 bit per field
                                      &(wahs[i]));
        (*wah_sizes)[i] = wahs_size[i];
        total_wah_size += wahs_size[i];
    }

    unsigned int W_i = 0;
    *W = (uint16_t *) malloc(total_wah_size*sizeof(uint16_t));
    for (i = 0; i < 4; i++) {
        for (j = 0; j < wahs_size[i]; j++) {
            (*W)[W_i] = wahs[i][j];
            W_i += 1;
        }
        free(wahs[i]);
    }

    return total_wah_size;
}
//}}}

//{{{ unsigned int ubin_to_bitmap_wah(unsigned int *U,
unsigned int ubin_to_bitmap_wah(unsigned int *U,
                                unsigned int U_len,
                                unsigned int num_fields,
                                unsigned int **W,
                                unsigned int **wah_sizes)
{
    unsigned int *B = NULL;
    // two bits per field
    unsigned int B_len = ubin_to_bitmap(U, U_len, num_fields*2, &B);
    unsigned int b_len = B_len / 4;  // size of each bitmap index

    *wah_sizes = (unsigned int *) malloc(4*sizeof(unsigned int));

    unsigned int *wahs[4],
                 wahs_size[4],
                 i,
                 j;

    unsigned int total_wah_size = 0;

    for (i = 0; i < 4; i++) {
        wahs_size[i] = ints_to_wah( (B + (i*b_len)), 
                                    b_len, 
                                    num_fields, // 1 bit per field
                                    &(wahs[i]));
        (*wah_sizes)[i] = wahs_size[i];
        total_wah_size += wahs_size[i];
    }

    free(B);

    unsigned int W_i = 0;
    *W = (unsigned int *) malloc(total_wah_size*sizeof(unsigned int));
    for (i = 0; i < 4; i++) {
        for (j = 0; j < wahs_size[i]; j++) {
            (*W)[W_i] = wahs[i][j];
            W_i += 1;
        }
        free(wahs[i]);
    }


    return total_wah_size;
}
//}}}

//{{{ unsigned int gt_records_ubin(struct ubin_file uf,
unsigned int gt_records_ubin(struct ubin_file uf,
                             unsigned int *record_ids,
                             unsigned int num_r,
                             unsigned int test_value,
                             unsigned int **R)
{
    return range_records_ubin(uf,
                              record_ids,
                              num_r,
                              test_value+1,
                              4,
                              R);

#if 0
    unsigned int num_output_ints = 1 + ((uf.num_fields - 1) / 32);

    unsigned int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    *R = (unsigned int *) malloc(num_output_ints*sizeof(unsigned int));
    unsigned int i,j;
    for (i = 0; i < num_output_ints; ++i)
        (*R)[i] = -1;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    unsigned int R_int_i, R_bit_i;
    unsigned int c_int_i, c_two_bit_i;

    for (i = 0; i < num_r; ++i) {
        // skip to the target record and read in the full record
        fseek(uf.file, uf.header_offset + // skip the record & field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

        R_int_i = 0;
        R_bit_i = 0;
        c_two_bit_i = 0;
        c_int_i = 0;

        for (j = 0; j < uf.num_fields; ++j) {
            // clear the bit

            if  ( !( test_value < ((c[c_int_i] >> (30 - c_two_bit_i*2)) & 3) ) )
                (*R)[R_int_i] = (*R)[R_int_i] & ~(1 << (31 - R_bit_i));

            R_bit_i += 1;
            if (R_bit_i == 32) {
                R_int_i += 1;
                R_bit_i = 0;
            }

            c_two_bit_i += 1;
            if (c_two_bit_i == 16) {
                c_two_bit_i = 0;
                c_int_i += 1;
            }
        }
    }

    free(c);

    return num_ints_per_record;
#endif
}
//}}}

//{{{ unsigned int print_ubin(struct ubin_file uf,
unsigned int print_ubin(struct ubin_file uf,
                        unsigned int *record_ids,
                        unsigned int num_r,
                        unsigned int format)
{
    unsigned int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    unsigned int i, j, k, num_printed = 0;

    if ( (num_r == 0) || (record_ids == NULL) ) {
        while (fread(c,sizeof(unsigned int),num_ints_per_record,uf.file)) {
            unsigned int printed_bits = 0;

            for (j = 0; j < num_ints_per_record; ++j) {
                if (j !=0)
                    printf(" ");

                if (format == 1) {
                    printf("%u", c[j]);
                } else if (format == 0) {
                    for (k = 0; k < 16; ++k) {
                        unsigned int bit = (c[j] >> (30 - 2*k)) & 3;
                        if (k !=0)
                            printf(" ");
                        printf("%u", bit);
                        printed_bits += 1;
                        if (printed_bits == uf.num_fields)
                            break;
                    }
                }
                num_printed += 1;
            }
            printf("\n");
        }
    } else {
        for (i = 0; i < num_r; ++i) {
            fseek(uf.file, uf.header_offset + 
                        record_ids[i]*num_bytes_per_record, 
                        SEEK_SET);
            fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

            unsigned int printed_bits = 0;

            for (j = 0; j < num_ints_per_record; ++j) {
                if (j !=0)
                    printf(" ");

                if (format == 1) {
                    printf("%u", c[j]);
                } else if (format == 0) {
                    for (k = 0; k < 16; ++k) {
                        unsigned int bit = (c[j] >> (30 - 2*k)) & 3;
                        if (k !=0)
                            printf(" ");
                        printf("%u", bit);
                        printed_bits += 1;
                        if (printed_bits == uf.num_fields)
                            break;
                    }
                }
                num_printed += 1;
            }
            printf("\n");
        }
    }

    free(c);

    return num_printed;
}
//}}}

//{{{ unsigned int print_by_name_ubin(char *ubin_file_name,
unsigned int print_by_name_ubin(char *ubin_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format)
{
    struct ubin_file uf = init_ubin_file(ubin_file_name);
    return print_ubin(uf, record_ids, num_r, format);
}
//}}} 

//{{{ unsigned int range_records_ubin(struct ubin_file uf,
unsigned int range_records_ubin(struct ubin_file uf,
                                unsigned int *record_ids,
                                unsigned int num_r,
                                unsigned int start_test_value,
                                unsigned int end_test_value,
                                unsigned int **R)
{
    unsigned int num_output_ints = 1 + ((uf.num_fields - 1) / 32);

    unsigned int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    *R = (unsigned int *) malloc(num_output_ints*sizeof(unsigned int));
    unsigned int i,j;
    for (i = 0; i < num_output_ints; ++i)
        (*R)[i] = -1;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    unsigned int R_int_i, R_bit_i;
    unsigned int c_int_i, c_two_bit_i;

    for (i = 0; i < num_r; ++i) {
        // skip to the target record and read in the full record
        fseek(uf.file, uf.header_offset + // skip the record & field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

        R_int_i = 0;
        R_bit_i = 0;
        c_two_bit_i = 0;
        c_int_i = 0;

        for (j = 0; j < uf.num_fields; ++j) {
            // clear the bit
            unsigned int val = (c[c_int_i] >> (30 - c_two_bit_i*2)) & 3;

            if (!(val >= start_test_value && val < end_test_value))
                (*R)[R_int_i] = (*R)[R_int_i] & ~(1 << (31 - R_bit_i));

            R_bit_i += 1;
            if (R_bit_i == 32) {
                R_int_i += 1;
                R_bit_i = 0;
            }

            c_two_bit_i += 1;
            if (c_two_bit_i == 16) {
                c_two_bit_i = 0;
                c_int_i += 1;
            }
        }
    }

    free(c);

    return num_ints_per_record;
}
//}}}

//{{{ unsigned int count_range_records_plt(struct plt_file pf,
unsigned int count_range_records_ubin(struct ubin_file uf,
                                      unsigned int *record_ids,
                                      unsigned int num_r,
                                      unsigned int start_test_value,
                                      unsigned int end_test_value,
                                      unsigned int **R)
{
    *R = (unsigned int *) calloc(uf.num_fields,sizeof(unsigned int));

    unsigned int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    unsigned int R_int_i, R_bit_i;
    unsigned int c_int_i, c_two_bit_i;

    unsigned int i,j;
    for (i = 0; i < num_r; ++i) {
        // skip to the target record and read in the full record
        fseek(uf.file, uf.header_offset + // skip the record & field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

        R_int_i = 0;
        R_bit_i = 0;
        c_two_bit_i = 0;
        c_int_i = 0;

        for (j = 0; j < uf.num_fields; ++j) {
            // clear the bit
            unsigned int val = (c[c_int_i] >> (30 - c_two_bit_i*2)) & 3;

            if ((val > start_test_value) && (val < end_test_value))
                (*R)[j] += 1;

            c_two_bit_i += 1;
            if (c_two_bit_i == 16) {
                c_two_bit_i = 0;
                c_int_i += 1;
            }
        }
    }

    free(c);

    return num_ints_per_record;
}
//}}}

//{{{ unsigned int gt_count_records_ubin(struct ubin_file uf,
unsigned int gt_count_records_ubin(struct ubin_file uf,
                                   unsigned int *record_ids,
                                   unsigned int num_r,
                                   unsigned int test_value,
                                   unsigned int **R)
{
    return count_range_records_ubin(uf,
                                   record_ids,
                                   num_r,
                                   test_value,
                                   4,
                                   R);
}
//}}}

//{{{ unsigned int print_ubin(struct ubin_file uf,
unsigned int convert_file_by_name_ubin_to_plt(char *ubin_in, char *plt_out)
{
    struct ubin_file uf = init_ubin_file(ubin_in);

    FILE *pf = fopen(plt_out,"w");

    if (!pf) {
        printf("Unable to open %s\n", plt_out);
        return 1;
    }

    fprintf(pf,"%u\n", uf.num_fields);
    fprintf(pf,"%u\n", uf.num_records);


    unsigned int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    unsigned int i, j, k, num_printed = 0;

    while (fread(c,sizeof(unsigned int),num_ints_per_record,uf.file)) {
        unsigned int printed_bits = 0;

        for (j = 0; j < num_ints_per_record; ++j) {
            if (j !=0)
                fprintf(pf," ");

            for (k = 0; k < 16; ++k) {
                unsigned int bit = (c[j] >> (30 - 2*k)) & 3;
                if (k !=0)
                    fprintf(pf," ");
                fprintf(pf,"%u", bit);
                printed_bits += 1;
                if (printed_bits == uf.num_fields)
                    break;
            }
            
           num_printed += 1;
        }
        fprintf(pf,"\n");
    }

    free(c);
    fclose(pf);
    fclose(uf.file);

    return num_printed;
}
//}}}
