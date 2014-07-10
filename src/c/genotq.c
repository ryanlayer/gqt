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
#include <math.h>
#include "genotq.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


static int test_1 = 0;
static char *line;
static int line_len;
static FILE *file;
static int read_state = 0;

//{{{ unsigned int bin_char_to_int(char *bin)
unsigned int bin_char_to_int(char *bin)
{
    unsigned int i = 0;
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

//{{{ struct plt_file init_plt_file(char *file_name)
struct plt_file init_plt_file(char *file_name)
{
    struct plt_file plt_f;

    plt_f.file = fopen(file_name, "r");

    if (!plt_f.file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        abort();
    }

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    int gt;

    read = getline(&line, &len, plt_f.file);

    plt_f.num_fields = atoi(line);


    read = getline(&line, &len, plt_f.file);

    plt_f.num_records = atoi(line);

    plt_f.header_offset = ftell(plt_f.file);

    free(line);

    return plt_f;
}
///}}}

//{{{ int or_records_plt(struct plt_file pf,
int or_records_plt(struct plt_file pf,
                   int *record_ids,
                   int num_r,
                   int *G)
{

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;

    long line_len = pf.num_fields*2*sizeof(char);
    int i,j;
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        for (j = 0; j < pf.num_fields; ++j)
            G[j] = G[j] | ((int)line[j*2] - 48);
    }

    free(line);

    return pf.num_fields;
}
//}}}

//{{{ int or_fields_plt(struct plt_file pf,
int or_fields_plt(struct plt_file pf,
                  int *field_ids,
                  int num_f,
                  int *G)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    long line_len = pf.num_fields*2*sizeof(char);
    int i,j;

    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);

    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        G[i] = 0;
        for (j = 0; j < num_f; ++j) 
            G[i] = G[i] | ((int)line[j*2] - 48);
    }

    free(line);

    return pf.num_records;
}
//}}}

//{{{unsigned int pack_2_bit_ints(int *ints, int num_ints)
unsigned int pack_2_bit_ints(int *ints, int num_ints)
{
    int i;
    unsigned int r = 0;
    for (i = 0; i < num_ints; ++i) {
       r = r | ((ints[i] & 3) << (30 - i*2));
    }
    
    return r;
}
//}}}

//{{{unsigned int *unpack_1_bit_ints(int packed_int)
int *unpack_1_bit_ints(unsigned int packed_ints)
{
    int *r = (int *) malloc (32*sizeof(int));

    int i;
    for (i = 0; i < 32; ++i) 
        r[i] = (packed_ints >> (31 - i)) & 1;
    
    return r;
}
//}}}

//{{{unsigned int *unpack_2_bit_ints(int packed_int)
int *unpack_2_bit_ints(unsigned int packed_ints)
{
    int *r = (int *) malloc (16*sizeof(int));

    int i;
    for (i = 0; i < 16; ++i) 
        r[i] = (packed_ints >> (30 - i*2)) & 3;
    
    return r;
}
//}}}

//{{{int convert_file_by_name_vcf_to_plt(char *in_file_name, char
int convert_file_by_name_vcf_to_plt(char *in_file_name,
                                    unsigned int num_fields,
                                    unsigned int num_records,
                                    char *out_file_name)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *pch;
    int gt;

    FILE *i_file = fopen(in_file_name, "r");
    if (!i_file) {
        fprintf(stderr, "Unable to open %s\n",in_file_name);
        return 1;
    }

    FILE *o_file = fopen(out_file_name, "w");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }


    fprintf(o_file, "%u\n%u\n", num_fields, num_records);

    while ((read = getline(&line, &len, i_file)) != -1) {
        if (line[0] != '#') {

            // Skip the first 9 fields
            
            int i;
            pch = strtok(line, "\t");
            for  (i =0; i < 8; ++i) {
                pch = strtok(NULL, "\t");
            }
            
            // get the first genotype
            pch = strtok(NULL, "\t");
            i = 0;
            while (pch != NULL) {
                if ((pch[0] == '0') && (pch[2] == '0'))
                    gt = 0;
                else if ((pch[0] == '1') && (pch[2] == '0'))
                    gt = 1;
                else if ((pch[0] == '0') && (pch[2] == '1'))
                    gt = 1;
                else if ((pch[0] == '1') && (pch[2] == '1'))
                    gt = 2;
                else
                    gt = 3;

                if (i != 0)
                    fprintf(o_file, " ");

                fprintf(o_file, "%d", gt);
                pch = strtok(NULL, "\t");
                i+=1;
            }
            fprintf(o_file, "\n");
        }
    }

    free(line);
    fclose(i_file);
    fclose(o_file);
    return 0;
}
//}}}

//{{{ void convert_file_by_name_plt_to_ubin(char *in_file_name, 
int convert_file_by_name_plt_to_ubin(char *in_file_name, char *out_file_name)
{
    struct plt_file pf = init_plt_file(in_file_name);
    int r = convert_file_plt_to_ubin(pf, out_file_name);
    fclose(pf.file);
    return r;
}

int convert_file_plt_to_ubin(struct plt_file pf, char *out_file_name)
{

    char *line = NULL;
    unsigned int i,j,count;
    size_t len;
    ssize_t read;
    int g[16];
    unsigned int *packed_ints;

    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First value is the number of fields per record
    fwrite(&(pf.num_fields), sizeof(int), 1, o_file);

    // Second value is the number of records
    fwrite(&(pf.num_records), sizeof(int), 1, o_file);


    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);
    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        //fprintf(stderr, "read:%lu\n", read);
        //fprintf(stderr, "%s", line);
        unsigned int num_packed_ints  = plt_line_to_packed_ints(line,
                                                                pf.num_fields,
                                                                &packed_ints);

        for (j = 0; j < num_packed_ints; ++j)
            fwrite(&(packed_ints[j]), sizeof(unsigned int), 1, o_file);

        free(packed_ints);
    }

    fclose(o_file);

    return 0;
}
//}}}

//{{{ unsigned int plt_line_to_packed_ints(char *line,
unsigned int plt_line_to_packed_ints(char *line,
                                     unsigned int num_fields, 
                                     unsigned int **packed_ints)
{
    unsigned int len = 1 + ((num_fields - 1) / 16);
    *packed_ints = (unsigned *) malloc((len)*sizeof(unsigned int));
    int i, two_bit_count, pack_int_count = 0;

    //printf("%s\n", line);
    
    int g[16];

    memset(g,0,sizeof(g));
    two_bit_count = 0;
    for (i = 0; i < num_fields; ++i) {
        g[i%16] = ((int)line[i*2] - 48);
        //printf("%d %d %c\n", i, i*2, line[i*2]);

        two_bit_count +=1;

        if (two_bit_count == 16) {
            int j;
            /*
            for (j = 0; j < two_bit_count; ++j) 
                printf("%d ", g[j]);
            printf("\n");
            */
            (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);
            pack_int_count += 1;
            memset(g,0,sizeof(g));
            two_bit_count = 0;
        }
        
    }


    /*
    for (i = 0; i < two_bit_count; ++i) 
        printf("%d ", g[i]);
    printf("\n");
    */

    if (two_bit_count > 0)
        (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);

    return len;
}
//}}}

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

//{{{ int or_records_ubin(struct ubin_file u_file, 
int or_records_ubin(struct ubin_file uf, 
                    int *record_ids,
                    int num_r,
                    unsigned int **G)
{

    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    // create enough space in r, and set all the values to zero
    *G = (unsigned int *) calloc(num_ints_per_record, sizeof(unsigned int));

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    int i,j;
    for (i = 0; i < num_r; ++i) {
        // skip to the target record and read in the full record
        fseek(uf.file, uf.header_offset + // skip the record & field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

        for (j = 0; j < num_ints_per_record; ++j)
            (*G)[j] = (*G)[j] | c[j];
    }

    free(c);

    return num_ints_per_record;
}
//}}}

//{{{ int or_fields_ubin(struct ubin_file u_file, 
/*
 *  This one is a little trickey.  To keep in like with the other opperations
 *  on ubin, we want G to be an array of packed ints, but since we are looking
 *  at each row and filling each two-bit value of the result we need to track both
 *  the current packed int and the possition within that int.
 *  The total number of two-bit values will be equal to the number of records.
 */
int or_fields_ubin(struct ubin_file uf, 
                   int *field_ids,
                   int num_f,
                   unsigned int **G)
{
    int num_ints_per_record = 1 + ((uf.num_fields - 1) / 16);

    int num_ints_all_recores = 1 + ((uf.num_records - 1) / 16);
    // create enough space in r, and set all the values to zero
    *G = (unsigned int *) calloc(num_ints_all_recores, sizeof(unsigned int));

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    // skip the record and field size field
    fseek(uf.file, uf.header_offset, SEEK_SET);

    int int_count = 0, two_bit_count = 0;
    unsigned int g, target_int;
    int field_int_id, field_offset;
    int i,j;
    for (i = 0; i < uf.num_records; ++i) {
        fread(c,sizeof(unsigned int),num_ints_per_record,uf.file);

        g = 0;
        for (j = 0; j < num_f; ++j) {
            field_int_id = field_ids[j]/16;
            field_offset = field_ids[j] - field_int_id*16;
            target_int = c[field_int_id];
            g = g | ((target_int >> (30 - field_offset*2)) & 3);
        }

        (*G)[int_count] = (*G)[int_count] | g << (30 - two_bit_count*2);

        two_bit_count += 1;

        if (two_bit_count == 16) {
            int_count += 1;
            two_bit_count = 0;
        }
    }

    if (two_bit_count > 0) {
        (*G)[int_count] = (*G)[int_count] | g << (30 - two_bit_count*2);
    }


    free(c);

    return num_ints_per_record;
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

//{{{unsigned int ints_to_rle(unsigned int *I, int I_len, unsigned int **O)
unsigned int ints_to_rle(unsigned int *I, int I_len, unsigned int **O)
{
    struct uint_ll *head=NULL,*tail=NULL;

    int i,
        j,
        curr_bit,
        ll_len = 0,
        last_bit = -1;
    unsigned int rle_v;
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

    *O = (unsigned int *) malloc(ll_len*sizeof(unsigned int));
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

//{{{unsinged int map_from_32_bits_to_15_bits(unsigned int *I,
unsigned int map_from_32_bits_to_15_bits(unsigned int *I,
                                         int I_len,
                                         unsigned int used_bits,
                                         uint16_t **O)
{
    unsigned int int_i, bit_i, group_i, in_group_i;
    unsigned int O_len =  (used_bits + 15 - 1)/ 15;
    //unsigned int O_len =  (I_len*32 + 31 - 1)/ 31;

    *O = (uint16_t *) calloc(O_len, sizeof(uint16_t));


    bit_i = 1;
    group_i = 0;
    in_group_i = 0;

    for (int_i = 0; int_i < I_len; ++int_i) {
        for ( ;bit_i<=32*(int_i+1); ++bit_i) {
            unsigned int bit = (I[int_i] >> (32 - (bit_i%32))) & 1;
            (*O)[group_i] = (*O)[group_i] + (bit << (14 - in_group_i));

            in_group_i += 1;

            if (bit_i % 15 == 0) {
                in_group_i = 0;
                group_i += 1;
            }
            if (bit_i == used_bits)
                break;
        }
        if (bit_i == used_bits)
            break;
    }
    return O_len;
}
//}}}

//{{{unsinged int map_from_32_bits_to_31_bits(unsigned int *I,
/* 
 * Take a list of 32 bit number and gives the list of 31-bit groups
 * (represented by 32-bit with left padding) ints 
 * Returns the number of 31-bit groups in O
 */
unsigned int map_from_32_bits_to_31_bits(unsigned int *I,
                                         int I_len,
                                         unsigned int used_bits,
                                         unsigned int **O)
{
    unsigned int int_i, bit_i, group_i, in_group_i;
    unsigned int O_len =  (used_bits + 31 - 1)/ 31;
    //unsigned int O_len =  (I_len*32 + 31 - 1)/ 31;

    *O = (unsigned int *) calloc(O_len, sizeof(unsigned int));


    bit_i = 1;
    group_i = 0;
    in_group_i = 0;

    for (int_i = 0; int_i < I_len; ++int_i) {
        for ( ;bit_i<=32*(int_i+1); ++bit_i) {
            unsigned int bit = (I[int_i] >> (32 - (bit_i%32))) & 1;
            (*O)[group_i] = (*O)[group_i] + (bit << (30 - in_group_i));

            in_group_i += 1;

            if (bit_i % 31 == 0) {
                in_group_i = 0;
                group_i += 1;
            }
            if (bit_i == used_bits)
                break;
        }
        if (bit_i == used_bits)
            break;
    }
    return O_len;
}
//}}}

//{{{ unsigned int wah_to_ints(unsigned int *W,
unsigned int wah_to_ints(unsigned int *W,
                         unsigned int W_len,
                         unsigned int **O)
{

    unsigned int wah_i;
    unsigned int num_bits = 0;

    for (wah_i = 0; wah_i < W_len; ++wah_i) {
        if (W[wah_i] >> 31 == 1) 
            num_bits += 31 * (W[wah_i] & 0x3fffffff); // zero out the fill bits
        else
            num_bits += 31;
    }

    unsigned int num_ints = (num_bits + 32 - 1) / 32;
    *O = (unsigned int *) malloc (num_ints * sizeof(unsigned int));


    unsigned int num_words,
                 word_i,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 int_i,
                 int_bit_i;

    int_bit_i = 1;
    int_i = 0;
    (*O)[int_i] = 0;
    for (wah_i = 0; wah_i < W_len; ++wah_i) {

        if (W[wah_i] >> 31 == 1) {
                num_words = (W[wah_i] & 0x3fffffff);
                fill_bit = (W[wah_i]>=0xC0000000?1:0);
                bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = W[wah_i];
        }

        for (word_i = 0; word_i < num_words; ++word_i) {
            for (bit_i = 0; bit_i < 31; ++bit_i) {
                bit = (bits >> (30 - bit_i)) & 1;
                //fprintf(stderr,"%u", bit);

                (*O)[int_i] += bit << (32 - int_bit_i);

                if (int_bit_i == 32) {
                    //fprintf(stderr,"\n%u\n",(*O)[int_i]);
                    int_i += 1;
                    int_bit_i = 0;
                    (*O)[int_i] = 0;
                }

                int_bit_i +=1;
            }
        }
    }
    
    return num_ints;
}
//}}}

//{{{ int append_bit_to_active_word(struct wah_active_word *a, int b)
int append_bit_to_active_word(struct wah_active_word *a, int b)
{
    // if a litteral, just add the next bit
    if (a->nbits < 31){ 
        a->value = (a->value<<1) + b;
        a->nbits += 1;
    // if a fill that isn't full, and value bit matches, add one
    //} else if ( (a->nbits > 31) && (((b<<30) & a->value) == b << 30)) {
    } else if ( (a->nbits > 31) && (((a->value >> 30) & 1) == b )) {
        if (a->nbits < pow(2,30) - 1) {
            a->nbits += 1;
            a->value += 1;
        } else {
            return 1; // can't add to current active
        }
    // if the litteral is full and they are all one value, make it a literal
    } else if ( (a->value == 0) || (a->value == pow(2,31) - 1)) {
        a->nbits = 32;
        a->value = ((2 + b) <<30) + 32;
    } else {
        return 1; // can't add to current active
    }

    return 0; // added to current active
}
//}}}

//{{{ int append_active_16word(struct wah_ll **A_head,
/*
 * This will return 1 if a new node is added to the active word linked list,
 * otherwise 0
 */
int append_active_16word(struct wah16_ll **A_head,
                         struct wah16_ll **A_tail,
                         struct wah16_active_word a)
{
    struct wah16_ll *n = (struct wah16_ll *)
        malloc(sizeof(struct wah16_ll));

    n->value = a;
    n->next = NULL;

    if (*A_head == NULL) {
        *A_head = n;
        *A_tail = n;
        return 1;
    } else if (a.value == 0) { // all zeros
        // The value on the tail is a litteral with all zeros
        if ( (*A_tail)->value.value == 0 ) {
            (*A_tail)->value.value = 0x8002;
            //(*A_tail)->value.value = 0x80000002;
            return 0;
        // The value on the tail is a fill of all zeros that is not full
        } else if ( ((*A_tail)->value.value >= 0x8000) &&
                    ((*A_tail)->value.value < 0xC000) ) {
            (*A_tail)->value.value += 1;
            return 0;
        } else { // the zeros cannot be added to the last active word
            (*A_tail)->next = n;
            *A_tail = n;
            return 1;
        }
    } else if (a.value == 0x7FFF) { // all ones
        if ( (*A_tail)->value.value == a.value ) {
            (*A_tail)->value.value = 0xC002;
            return 0;
        } else if ( (*A_tail)->value.value >= 0xC000) {
            (*A_tail)->value.value += 1;
            return 0;
        } else {
            (*A_tail)->next = n;
            *A_tail = n;
            return 1;
        }
    } else {
        (*A_tail)->next = n;
        *A_tail = n;
        return 1;
    }
}
//}}}

//{{{ int append_active_word(struct wah_ll **A_head,
/*
 * This will return 1 if a new node is added to the active word linked list,
 * otherwise 0
 */
int append_active_word(struct wah_ll **A_head,
                       struct wah_ll **A_tail,
                       struct wah_active_word a)
{
    struct wah_ll *n = (struct wah_ll *)
        malloc(sizeof(struct wah_ll));

    n->value = a;
    n->next = NULL;

    if (*A_head == NULL) {
        *A_head = n;
        *A_tail = n;
        return 1;
    } else if (a.value == 0) { // all zeros
        // The value on the tail is a litteral with all zeros
        if ( (*A_tail)->value.value == 0 ) {
            (*A_tail)->value.value = 0x80000002;
            //(*A_tail)->value.value = (1<<31) + 62;
            //(*A_tail)->value.nbits = 62;
            return 0;
        // The value on the tail is a fill of all zeros that is not full
        /*
        } else if (  ((*A_tail)->value.value >> 30 == 2) &&
                    (((*A_tail)->value.value << 2) >> 2 < (pow(2,30)-31)) ) {
        */
        } else if ( ((*A_tail)->value.value >= 0x80000000) &&
                    ((*A_tail)->value.value < 0xC0000000) ) {
            (*A_tail)->value.value += 1;
            //(*A_tail)->value.value += 31;
            //(*A_tail)->value.nbits += 31;
            return 0;
        } else { // the zeros cannot be added to the last active word
            (*A_tail)->next = n;
            *A_tail = n;
            return 1;
        }
    //} else if (a.value == (pow(2,31) - 1)) { // all ones
    } else if (a.value == 0x7FFFFFFF) { // all ones
        if ( (*A_tail)->value.value == a.value ) {
            (*A_tail)->value.value = 0xC0000002;
            //(*A_tail)->value.value = (3<<30) + 62;
            //(*A_tail)->value.nbits = 62;
            return 0;
        /*
        } else if (  ((*A_tail)->value.value >> 30 == 3) &&
                    (((*A_tail)->value.value << 2) >> 2 < (pow(2,30)-31)) ) {
        */
        } else if ( (*A_tail)->value.value >= 0xC0000000) {
            (*A_tail)->value.value += 1;
            //(*A_tail)->value.value += 31;
            //(*A_tail)->value.nbits += 31;
            return 0;
        } else {
            (*A_tail)->next = n;
            *A_tail = n;
            return 1;
        }
    } else {
        (*A_tail)->next = n;
        *A_tail = n;
        return 1;
    }
}
//}}}

//{{{ int append_fill_word(struct wah_ll **A_head,
/*
 * The fill_size is the number of words / 31*fill_size number of bits
 */
int append_fill_word(struct wah_ll **A_head,
                     struct wah_ll **A_tail,
                     int fill_bit,
                     unsigned int fill_size)
{

    // if it is a fill, the bit matches, and there is room 
    if ( (*A_head != NULL) &&
         ((*A_tail)->value.value >> 30 == 2 + fill_bit) &&//fill & bit match
         (((*A_tail)->value.value >> 2)+ fill_size < 0x3fffffff) ) {//room

        (*A_tail)->value.value += fill_size;
        return 0;

    } else {
        struct wah_ll *n = (struct wah_ll *)
                malloc(sizeof(struct wah_ll));

        n->next = NULL;

        if (fill_size > 1)
            n->value.value = ((2 + fill_bit) << 30) + fill_size;
        else { 
            n->value.nbits = 31;
            n->value.value = (fill_bit?0x7FFFFFFF:0);
        }

        if (*A_head == NULL)
            *A_head = n;
        else
            (*A_tail)->next = n;

        *A_tail = n;

        return 1;
    }
    
    return -1;
} 
//}}}

//{{{ int ints_to_wah16(unsigned int *I,
unsigned int ints_to_wah16(unsigned int *I,
                           int I_len,
                           unsigned int used_bits,
                           uint16_t **W)
{
    unsigned int W_len;
    uint16_t *O;
    // split the intput up int to 31-bit groups
    unsigned int O_len = map_from_32_bits_to_15_bits(I, I_len, used_bits, &O);

    // build the WAH list
    struct wah16_ll *A_head = NULL,
                    *A_tail = NULL,
                    *A_curr;
    int i,c = 0;
    struct wah16_active_word a;
    a.nbits=15;
    for (i = 0; i < O_len; ++i) {
        a.value = O[i];
        c += append_active_16word(&A_head,&A_tail,a);
    }

    free(O);

    // Move the linked list to an array
    W_len = c;
    *W = (uint16_t *) malloc(W_len * sizeof(uint16_t));

    A_curr = A_head;
    struct wah16_ll *A_tmp;
    i = 0;
    while (A_curr != NULL) {
        (*W)[i] = A_curr->value.value;
        i += 1;
        A_tmp = A_curr;
        A_curr = A_tmp->next;
        free(A_tmp);
    }

    return W_len;
}
//}}}

//{{{ int ints_to_wah(unsigned int *I,
unsigned int ints_to_wah(unsigned int *I,
                         int I_len,
                         unsigned int used_bits,
                         unsigned int **W)
{
    unsigned int W_len;
    unsigned int *O;
    // split the intput up int to 31-bit groups
    unsigned int O_len = map_from_32_bits_to_31_bits(I, I_len, used_bits, &O);

    // build the WAH list
    struct wah_ll *A_head = NULL,
                   *A_tail = NULL,
                   *A_curr;
    int i,c = 0;
    struct wah_active_word a;
    a.nbits=31;
    for (i = 0; i < O_len; ++i) {
        a.value = O[i];
        c += append_active_word(&A_head,&A_tail,a);
    }

    free(O);

    // Move the linked list to an array
    W_len = c;
    *W = (unsigned int *) malloc(W_len * sizeof(unsigned int));

    A_curr = A_head;
    struct wah_ll *A_tmp;
    i = 0;
    while (A_curr != NULL) {
        (*W)[i] = A_curr->value.value;
        i += 1;
        A_tmp = A_curr;
        A_curr = A_tmp->next;
        free(A_tmp);
    }

    return W_len;
}
//}}}

//{{{ struct wah_run init_wah_run(unsigned int *words,
struct wah_run init_wah_run(unsigned int *words,
                            unsigned int len){
    struct wah_run r;
    r.words = words;
    r.word_i = 0;
    r.len = len;
    r.fill = 0;
    r.num_words = 0;
    r.is_fill = 0;

    return r;
}
//}}}

//{{{ void wah_run_decode(struct wah_run *r)
void wah_run_decode(struct wah_run *r)
{
    if (r->words[r->word_i] > 0x7FFFFFFF) {
        r->fill = (r->words[r->word_i]>=0xC0000000?0x7FFFFFFF:0);
        r->num_words = r->words[r->word_i] & 0x3FFFFFFF;
        r->is_fill= 1;
        r->fill_bit = (r->words[r->word_i]>=0xC0000000?1:0);
    } else {
        r->num_words = 1;
        r->is_fill= 0;
    }
}
//}}}

//{{{ void wah_or(struct wah_run *x,
unsigned int  wah_or(struct wah_run *x,
                     struct wah_run *y,
                     unsigned int **O)
{
    struct wah_active_word a;
    unsigned int num_words;
    struct wah_ll *Z_head = NULL,
                  *Z_tail = NULL;
    int Z_len = 0;

    //xrun.it = x.vec.begin(); xrun.decode();
    wah_run_decode(x);
    //yrun.it = y.vec.begin(); yrun.decode();
    wah_run_decode(y);

    // WHILE (x.vec and y.vec are not exhausted) 
    while ( (x->word_i < x->len) && (y->word_i < y->len) ){
        // IF (xrun.nWords == 0) ++xrun.it, xrun.decode();
        if (x->num_words == 0) {
            x->word_i += 1;
            if (x->word_i < x->len)
                wah_run_decode(x);
        }

        // IF (yrun.nWords == 0) ++yrun.it, yrun.decode();
        if (y->num_words == 0) {
            y->word_i += 1;
            if (y->word_i < y->len)
                wah_run_decode(y);
        }

        if ( (x->word_i >= x->len) && (y->word_i >= y->len) )
            break;
        else if (x->word_i >= x->len) {
            fprintf(stderr, "X ended before Y\n");
            abort();
        } else if (y->word_i >= y->len) {
            fprintf(stderr, "Y ended before X\n");
            abort();
        }


        //   IF (xrun.isFill)
        if (x->is_fill == 1) {
            // IF (yrun.isFill)
            if (y->is_fill == 1) {
                //fprintf(stderr,"X Fill\tY Fill\n");
                //nWords = min(xrun.nWords, yrun.nWords)
                //z.appendFill(nWords, (*(xrun.it) ◦ *(yrun.it))),
                //xrun.nWords -= nWords, yrun.nWords -= nWords;
                num_words = MIN(x->num_words, y->num_words);
                Z_len += append_fill_word(&Z_head,
                                          &Z_tail,
                                          (x->fill_bit | y->fill_bit),
                                          num_words);
                x->num_words -= num_words;
                y->num_words -= num_words;
            // ELSE
            } else {
                //fprintf(stderr,"X Fill\tY Litt\n");
                //z.active.value = xrun.fill ◦ *yrun.it
                //z.appendLiteral(),
                //-- xrun.nWords, yrun.nWords = 0;
                a.nbits = 31;
                a.value = y->words[y->word_i] | x->fill;
                Z_len += append_active_word(&Z_head,
                                            &Z_tail,
                                            a);
                x->num_words -= 1;
                y->num_words = 0;
            }
        //ELSEIF (yrun.isFill)
        } else if (y->is_fill == 1) {
            //fprintf(stderr,"X Litt\tY Fill\n");
            //z.active.value = yrun.fill ◦ *xrun.it,
            //z.appendLiteral(),
            //-- yrun.nWords, xrun.nWords = 0;
            a.nbits = 31;
            a.value = x->words[x->word_i] | y->fill;

            Z_len += append_active_word(&Z_head,
                                        &Z_tail,
                                        a);
            y->num_words -= 1;
            x->num_words = 0;
        //ELSE
        } else {
            //fprintf(stderr,"X Litt\tY Litt\n");
            a.nbits = 31;
            a.value = x->words[x->word_i] | y->words[y->word_i];
            Z_len += append_active_word(&Z_head,
                                        &Z_tail,
                                        a);
            y->num_words = 0;
            x->num_words = 0;
        }
    }

    *O = (unsigned int *) malloc(Z_len * sizeof(unsigned int));

    struct wah_ll *Z_tmp, *Z_curr = Z_head;
    int i = 0;
    while (Z_curr != NULL) {
        (*O)[i] = Z_curr->value.value;
        i += 1;
        Z_tmp = Z_curr;
        Z_curr = Z_tmp->next;
        free(Z_tmp);
    }

    return Z_len;
}
//}}}

//{{{ void wah_or(struct wah_run *x,
unsigned int  wah_in_place_or(unsigned int *r_wah,
                              unsigned int r_wah_size,
                              unsigned int *wah,
                              unsigned int wah_size)
{

    unsigned int r_wah_i = 0;
    unsigned int wah_i, fill_size, wah_v, end;
    for (wah_i = 0; wah_i < wah_size; ++wah_i)
    {
        wah_v = wah[wah_i];
        // is the current word a fill
        if (wah_v >= 0x80000000) {
            fill_size = wah_v & 0x3fffffff;
            if (wah_v >> 30 == 3) {
                // fill of 1s
                end =  r_wah_i + fill_size;
                for ( ; r_wah_i < end; ++r_wah_i)
                    r_wah[r_wah_i] = 0x7fffffff;
            } else {
                // fill of 0s
                r_wah_i += fill_size;
            }
        } else {
            r_wah[r_wah_i] = r_wah[r_wah_i] | wah[wah_i];
            r_wah_i += 1;
        }
    }

    return r_wah_size;
}
//}}}

//{{{ void wah_or(struct wah_run *x,
unsigned int  wah_compressed_in_place_or(unsigned int *r_wah,
                                         unsigned int r_wah_size,
                                         unsigned int *wah,
                                         unsigned int wah_size)
{
    unsigned int wah_i, wah_v, wah_fill_size, wah_fill_value,
                 r_wah_i, r_wah_v, r_wah_fill_size, r_wah_fill_value,
                 end, num_words;

    r_wah_i = 0;

    for (wah_i = 0; wah_i < wah_size; ++wah_i)
    {
        wah_v = wah[wah_i];
        r_wah_v = r_wah[r_wah_i];

        if (wah_v >= 0x80000000) {
            wah_fill_value = (wah_v >> 30) & 1;
            wah_fill_size = (wah_v & 0x3fffffff);

            while (wah_fill_size > 0) {

                if (r_wah_v >= 0x80000000) { // r_wah is a fill

                    r_wah_fill_value = (r_wah_v >> 30) & 1;
                    r_wah_fill_size = (r_wah_v & 0x3fffffff);

                    // make a new fill based on the smaller one
                    num_words = MIN(wah_fill_size, r_wah_fill_size);

                    r_wah[r_wah_i] = (1 << 31) + 
                                     ((r_wah_fill_value | 
                                        wah_fill_value) << 30) + 
                                     num_words;

                    r_wah_fill_size -= num_words;
                    wah_fill_size -= num_words;

                    // save any values left on the end of r_wah run
                    if (r_wah_fill_size > 0) {
                        if (r_wah_fill_size == 1) {
                            // we no longer have a fill, write a literal
                            if (r_wah_fill_value == 1) //all ones
                                r_wah[r_wah_i + num_words] = 0x7fffffff;
                            else  // all zeros
                                r_wah[r_wah_i + num_words] = 0;
                        } else { 
                            // we still have a fill, write it
                            r_wah[r_wah_i + num_words] = 
                                    (1 << 31) + 
                                    (r_wah_fill_value << 30) + 
                                    r_wah_fill_size; 
                        }
                    }

                    r_wah_i += num_words;
                } else { // r_wah is a literal
                    if (wah_fill_value == 1)
                        r_wah[r_wah_i] = 0x7fffffff;
                    r_wah_i += 1;
                }

                r_wah_v = r_wah[r_wah_i];
            }
        } else if ( r_wah_v >= 0x80000000) {
            r_wah_fill_value = (r_wah_v >> 30) & 1;

            //wah is not a fill and r_wah is a fill
            //update the current word in r_wah 
            if (r_wah_fill_value == 1) {
                //fill is a one
                r_wah[r_wah_i] = 0x7fffffff;
            } else {
                //fill is a zero
                r_wah[r_wah_i] = wah_v;
            }

            // we just took one word off the r_wah fill
            r_wah_fill_size = (r_wah_v & 0x3fffffff) - 1;

            if (r_wah_fill_size == 1) {
                // we no longer have a fill, write a literal
                if (r_wah_fill_value == 1) //all ones
                    r_wah[r_wah_i + 1] = 0x7fffffff;
                else  // all zeros
                    r_wah[r_wah_i + 1] = 0;
            } else { 
                // we still have a fill, write it
                r_wah[r_wah_i + 1] = (1 << 31) + 
                                     (r_wah_fill_value << 30) + 
                                     r_wah_fill_size; 
            }
            r_wah_i += 1;

        } else {
            r_wah[r_wah_i] = r_wah[r_wah_i] | wah[wah_i];
            r_wah_i += 1;
        }
    }

    return r_wah_size;
}
//}}}

//{{{ void wah_and(struct wah_run *x,
unsigned int  wah_in_place_and(unsigned int *r_wah,
                               unsigned int r_wah_size,
                               unsigned int *wah,
                               unsigned int wah_size)
{

    unsigned int r_wah_i = 0;
    unsigned int wah_i, fill_size, wah_v, end;
    for (wah_i = 0; wah_i < wah_size; ++wah_i)
    {
        wah_v = wah[wah_i];
        // is the current word a fill
        if (wah_v >= 0x80000000) {
            fill_size = wah_v & 0x3fffffff;
            if (wah_v >> 30 == 3) {
                // fill of 1s
                r_wah_i += fill_size;
            } else {
                // fill of 0s
                end =  r_wah_i + fill_size;
                for ( ; r_wah_i < end; ++r_wah_i)
                    r_wah[r_wah_i] = 0;
            }
        } else {
            r_wah[r_wah_i] = r_wah[r_wah_i] & wah[wah_i];
            r_wah_i += 1;
        }
    }

    return r_wah_i;
}
//}}}

//{{{ void wah_and(struct wah_run *x,
unsigned int  wah_and(struct wah_run *x,
                     struct wah_run *y,
                     unsigned int **O)
{
    struct wah_active_word a;
    unsigned int num_words;
    struct wah_ll *Z_head = NULL,
                  *Z_tail = NULL;
    int Z_len = 0;

    //xrun.it = x.vec.begin(); xrun.decode();
    wah_run_decode(x);
    //yrun.it = y.vec.begin(); yrun.decode();
    wah_run_decode(y);

    // WHILE (x.vec and y.vec are not exhausted) 
    while ( (x->word_i < x->len) && (y->word_i < y->len) ){
        // IF (xrun.nWords == 0) ++xrun.it, xrun.decode();
        if (x->num_words == 0) {
            x->word_i += 1;
            if (x->word_i < x->len)
                wah_run_decode(x);
        }

        // IF (yrun.nWords == 0) ++yrun.it, yrun.decode();
        if (y->num_words == 0) {
            y->word_i += 1;
            if (y->word_i < y->len)
                wah_run_decode(y);
        }

        if ( (x->word_i >= x->len) && (y->word_i >= y->len) )
            break;
        else if (x->word_i >= x->len) {
            fprintf(stderr, "X ended before Y\n");
            abort();
        } else if (y->word_i >= y->len) {
            fprintf(stderr, "Y ended before X\n");
            abort();
        }


        //   IF (xrun.isFill)
        if (x->is_fill == 1) {
            // IF (yrun.isFill)
            if (y->is_fill == 1) {
                //fprintf(stderr,"X Fill\tY Fill\n");
                //nWords = min(xrun.nWords, yrun.nWords)
                //z.appendFill(nWords, (*(xrun.it) ◦ *(yrun.it))),
                //xrun.nWords -= nWords, yrun.nWords -= nWords;
                num_words = MIN(x->num_words, y->num_words);
                Z_len += append_fill_word(&Z_head,
                                          &Z_tail,
                                          (x->fill_bit & y->fill_bit),
                                          num_words);
                x->num_words -= num_words;
                y->num_words -= num_words;
            // ELSE
            } else {
                //fprintf(stderr,"X Fill\tY Litt\n");
                //z.active.value = xrun.fill ◦ *yrun.it
                //z.appendLiteral(),
                //-- xrun.nWords, yrun.nWords = 0;
                a.nbits = 31;
                a.value = y->words[y->word_i] & x->fill;
                Z_len += append_active_word(&Z_head,
                                            &Z_tail,
                                            a);
                x->num_words -= 1;
                y->num_words = 0;
            }
        //ELSEIF (yrun.isFill)
        } else if (y->is_fill == 1) {
            //fprintf(stderr,"X Litt\tY Fill\n");
            //z.active.value = yrun.fill ◦ *xrun.it,
            //z.appendLiteral(),
            //-- yrun.nWords, xrun.nWords = 0;
            a.nbits = 31;
            a.value = x->words[x->word_i] & y->fill;

            Z_len += append_active_word(&Z_head,
                                        &Z_tail,
                                        a);
            y->num_words -= 1;
            x->num_words = 0;
        //ELSE
        } else {
            //fprintf(stderr,"X Litt\tY Litt\n");
            a.nbits = 31;
            a.value = x->words[x->word_i] & y->words[y->word_i];
            Z_len += append_active_word(&Z_head,
                                        &Z_tail,
                                        a);
            y->num_words = 0;
            x->num_words = 0;
        }
    }

    *O = (unsigned int *) malloc(Z_len * sizeof(unsigned int));

    struct wah_ll *Z_tmp, *Z_curr = Z_head;
    int i = 0;
    while (Z_curr != NULL) {
        (*O)[i] = Z_curr->value.value;
        i += 1;
        Z_tmp = Z_curr;
        Z_curr = Z_tmp->next;
        free(Z_tmp);
    }

    return Z_len;
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
    unsigned int *B;
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

//{{{ unsigned int plt_to_bitmap_wah(unsigned int *U,
unsigned int plt_to_bitmap_wah(char *plt,
                               unsigned int plt_len,
                               unsigned int **W,
                               unsigned int **wah_sizes)
{

    unsigned int *ubin;
    unsigned int ubin_len = plt_line_to_packed_ints(plt, plt_len, &ubin);

    unsigned int wah_len = ubin_to_bitmap_wah(ubin,
                                              ubin_len,
                                              plt_len, // one bit per field
                                              W,
                                              wah_sizes);
    free(ubin);

    return wah_len;
}
//}}}

//{{{ unsigned int convert_file_by_name_ubin_to_wahbm(char *ubin_in, 
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

//{{{ struct wah_file init_wah_file(char *file_name)
struct wah_file init_wah_file(char *file_name)
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
            malloc(sizeof (unsigned int)*wf.num_records);

    unsigned int i;
    for (i = 0; i < wf.num_records; ++i)
        fread(&(wf.record_offsets[i]),sizeof(unsigned int),1,wf.file);


    wf.header_offset = ftell(wf.file);

    return wf;
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

//{{{ unsigned int get_wah_bitmap(struct wah_file wf,
unsigned int get_wah_bitmap_in_place(struct wah_file wf,
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


    //*wah_bitmap = (unsigned int *) malloc(sizeof(unsigned int)*wah_size);
    fseek(wf.file, wah_offset, SEEK_SET);
    fread(*wah_bitmap,sizeof(unsigned int),wah_size,wf.file);

    return wah_size;
}
//}}}

//{{{ unsigned int get_wah_record(struct wah_file wf,
unsigned int get_wah_record(struct wah_file wf,
                            unsigned int wah_record,
                            unsigned int **wah)
{
    // get the size of the WAH-encoded bitmap
    unsigned int wah_size = 0, wah_offset = 0;
    if ( wah_record == 0) {
        wah_size = wf.record_offsets[wah_record];
        wah_offset = wf.header_offset;
    } else {
        wah_size = wf.record_offsets[wah_record] - 
                   wf.record_offsets[wah_record - 1];

        wah_offset = wf.header_offset +
                     sizeof(unsigned int) * 
                        (wf.record_offsets[wah_record] - wah_size);
    }

    *wah = (unsigned int *) malloc(sizeof(unsigned int)*wah_size);
    fseek(wf.file, wah_offset, SEEK_SET);
    fread(*wah,sizeof(unsigned int),wah_size,wf.file);

    return wah_size;
}
//}}}

//{{{ unsigned int get_plt_record(struct plt_file pf,
unsigned int get_plt_record(struct plt_file pf,
                            unsigned int plt_record,
                            unsigned int **plt)
{

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *plt = (unsigned int *) calloc(num_ints_per_record, sizeof(unsigned int));

    long line_len = pf.num_fields*2*sizeof(char);

    fseek(pf.file, pf.header_offset + line_len*plt_record, SEEK_SET);
    read = getline(&line, &len, pf.file);

    unsigned int plt_size = plt_line_to_packed_ints(line,
                                                    pf.num_fields, 
                                                    plt);
    free(line);
    return plt_size;
}
//}}}

//{{{ unsigned int range_records_plt(struct plt_file pf,
unsigned int range_records_plt(struct plt_file pf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int **R)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *R = (unsigned int *) malloc(num_ints_per_record*sizeof(unsigned int));

    for (i = 0; i < num_ints_per_record; ++i)
        (*R)[i] = -1;

    unsigned int int_i, bit_i;


    long line_len = pf.num_fields*2*sizeof(char);
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        int_i = 0;
        bit_i = 0;
        for (j = 0; j < pf.num_fields; ++j) {
            // clear the bit if fails criteria
            unsigned int val = (unsigned int)line[j*2] - 48;            
            if (!(val >= start_test_value && val < end_test_value))
                (*R)[int_i] = (*R)[int_i] & ~(1 << (31 - bit_i));

            bit_i += 1;
            if (bit_i == 32) {
                int_i += 1;
                bit_i = 0;
            }
        }
    }

    free(line);

    return num_ints_per_record;
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int eq_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value, test_value+1, R);
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int ne_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *R = (unsigned int *) malloc(num_ints_per_record*sizeof(unsigned int));

    for (i = 0; i < num_ints_per_record; ++i)
        (*R)[i] = -1;

    unsigned int int_i, bit_i;


    long line_len = pf.num_fields*2*sizeof(char);
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        int_i = 0;
        bit_i = 0;
        for (j = 0; j < pf.num_fields; ++j) {
            // clear the bit
            unsigned int val = ((unsigned int)line[j*2] - 48);
            if  (!(val != test_value))
                (*R)[int_i] = (*R)[int_i] & ~(1 << (31 - bit_i));

            bit_i += 1;
            if (bit_i == 32) {
                int_i += 1;
                bit_i = 0;
            }
        }
    }

    free(line);

    return num_ints_per_record;
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int gt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value+1, 4, R);
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int gte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value, 4, R);
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int lt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, 0, test_value, R);
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int lte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, 0, test_value+1, R);
}
//}}}

//{{{ unsigned int gt_records_ubin(struct ubin_file uf,
unsigned int gt_records_ubin(struct ubin_file uf,
                             unsigned int *record_ids,
                             unsigned int num_r,
                             unsigned int test_value,
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

//{{{ unsigned int range_records_in_place_wahbm(struct wah_file wf,
unsigned int range_records_in_place_wahbm(struct wah_file wf,
                                          unsigned int *record_ids,
                                          unsigned int num_r,
                                          unsigned int start_test_value,
                                          unsigned int end_test_value,
                                          unsigned int **R) 

{

    unsigned int max_wah_size = (wf.num_fields + 31 - 1)/ 31;
    unsigned int *record_new_bm = (unsigned int *)
                        malloc(sizeof(unsigned int)*max_wah_size);

    unsigned int *or_result_bm = (unsigned int *)
                        malloc(sizeof(unsigned int)*max_wah_size);
    unsigned int *and_result_bm = (unsigned int *)
                        malloc(sizeof(unsigned int)*max_wah_size);
    unsigned int and_result_bm_size, record_new_bm_size, or_result_bm_size;
    unsigned int i,j;
    for (i = 0; i < max_wah_size; ++i)
        and_result_bm[i] = 0x7fffffff;

    for (i = 0; i < num_r; ++i) {
        // or the appropriate bitmaps
        memset(or_result_bm, 0, sizeof(unsigned int)*max_wah_size);

        for (j = start_test_value; j < end_test_value; ++j) {

            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);

            or_result_bm_size = wah_in_place_or(or_result_bm,
                                                max_wah_size,
                                                record_new_bm,
                                                record_new_bm_size); 
        }

        // and 
        and_result_bm_size = wah_in_place_and(and_result_bm,
                                              max_wah_size,
                                              or_result_bm,
                                              or_result_bm_size);

    }

    free(record_new_bm);
    free(or_result_bm);

    *R = and_result_bm;
    return and_result_bm_size;
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

//{{{ unsigned int eq_records_wahbm(struct wah_file wf,
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

//{{{ unsigned int gt_records_in_place_wahbm(struct wah_file wf,
unsigned int gt_records_in_place_wahbm(struct wah_file wf,
                                       unsigned int *record_ids,
                                       unsigned int num_r,
                                       unsigned int test_value,
                                       unsigned int **R) 

{
    // TODO: need constants for upper bound.
    return range_records_in_place_wahbm(wf,
                                        record_ids,
                                        num_r,
                                        test_value+1,
                                        4,
                                        R);
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

//{{{ unsigned int lt_records_wahbm(struct wah_file wf,
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

//{{{ unsigned int print_plt(struct plt_file pf,
unsigned int print_plt(struct plt_file pf,
                       unsigned int *record_ids,
                       unsigned int num_r)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;
    unsigned int num_printed = 0;


    long line_len = pf.num_fields*2*sizeof(char);

    if ( (num_r == 0) || record_ids == NULL ) {
        while ((read = getline(&line, &len, pf.file)) != -1) {
            printf("%s", line);
            num_printed += 1;
        }
    } else {
        for (i = 0; i < num_r; ++i) {
            fseek(pf.file,
                  pf.header_offset + line_len*record_ids[i],
                  SEEK_SET);

            read = getline(&line, &len, pf.file);
            printf("%s", line);
            num_printed += 1;
        }
    }

    free(line);

    return num_printed;
}
//}}}

//{{{unsigned int print_by_name_plt(char *pf_file_name,
unsigned int print_by_name_plt(char *pf_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r)
{
    struct plt_file pf = init_plt_file(pf_file_name);
    return print_plt(pf, record_ids, num_r); 

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

//{{{ unsigned int print_wahbm(struct ubin_file wf,
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

//{{{ unsigned int print_wah(struct ubin_file wf,
unsigned int print_wah(struct wah_file wf,
                       unsigned int *record_ids,
                       unsigned int num_r,
                       unsigned int format)
{
    unsigned int i,j,k,wah_size,printed_bits,to_print = num_r;
    unsigned int *wah = NULL;

    if (num_r == 0)
        to_print = wf.num_records;

    for (i = 0; i < to_print; ++i) {

        // get the compressed bitmap
        if (num_r > 0)
            wah_size = get_wah_record(wf,
                                     record_ids[i],
                                     &wah);
        else
            wah_size = get_wah_record(wf,
                                     i,
                                     &wah);

        // decompress 
        unsigned int *ints = NULL;
        unsigned int ints_size = wah_to_ints(wah,wah_size,&ints);
        printed_bits = 0;

        for (j = 0; j < ints_size; ++j) {
            if (j !=0)
                printf(" ");
            for (k = 0; k < 16; ++k) {
                unsigned int val = (ints[j] >> (30 - 2*k)) & 3;
                if (k !=0)
                    printf(" ");
                printf("%u", val);
                printed_bits += 1;
                if (printed_bits == wf.num_fields)
                    break;
            }
            if (printed_bits == wf.num_fields)
                break;
        }
        printf("\n");

        free(ints);
        ints = NULL;
        free(wah);
        wah = NULL;
    }

    return to_print;
}
//}}}

//{{{ unsigned int print_by_name_wah(char *wahbm_file_name,
unsigned int print_by_name_wah(char *wahbm_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r,
                               unsigned int format)
{
    struct wah_file wf = init_wah_file(wahbm_file_name);
    return print_wah(wf, record_ids, num_r, format);
}
//}}} 

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//{{{ void init_int_genotype_reader(char *file_name, int num_gt)
void init_int_genotype_reader(char *file_name, int num_gt)
{

    line_len = num_gt*2 + 1;
    line = (char *) malloc(line_len*sizeof(char));
    file = fopen(file_name, "r");
    read_state = 0;
}
//}}}

//{{{void destroy_int_genotype_reader()
void destroy_int_genotype_reader()
{
    fclose(file);
    free(line);
}
//}}}

//{{{ int get_next_int_genotype(int *line_num, int *gt_num, int *gt)
int get_next_int_genotype(int *line_num, int *gt_num, int *gt)
{
    // read state: 0 get a line
    // read state: 1 get first char
    // read state: 2 get next char
    // read state: 3 EOF

    char *pch;

    if (read_state == 0) {
        char *r = fgets(line, line_len, file);
        if (r == NULL) {
            read_state = 3;
            return 1;
        } else {
            read_state = 1;
        }
    }

    if (read_state == 1) { 
        pch = strtok(line, " ");
       *gt_num = 0;

        read_state = 2;
    } else if (read_state == 2) {
        pch = strtok(NULL, " ");
       *gt_num = *gt_num + 1;
    }

    if (pch != NULL) {
        *gt = atoi(pch);
        return 0;
    } else {
        read_state = 0;
        *line_num = *line_num + 1;
        return get_next_int_genotype(line_num, gt_num, gt);
    }
}
//}}}

//{{{struct uint_file init_uint_file(char *file_name, 
struct uint_file init_uint_file(char *file_name, 
                                int num_records,
                                int num_fields)
{
    struct uint_file uf;
    uf.line_len = num_fields*2+1;
    uf.file = fopen(file_name, "r");

    if (uf.file == NULL) {
        printf("Could not open %s\n", file_name);
        abort();
    }

    uf.num_records = num_records;
    uf.num_fields = num_fields;

    return uf;
}
//}}}

//{{{int or_uint_records(struct uint_file u_file, 
int or_uint_records(struct uint_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r)
{
    *r = (unsigned int *) calloc(u_file.num_fields, sizeof(unsigned int));

    int a = 0;
    int i,j;
    char *pch, *ret;
    char *l = (char *) malloc(u_file.line_len*sizeof(char));
    for (i = 0; i < num_r; ++i) {
        fseek(u_file.file,u_file.num_fields*2 * record_ids[i],SEEK_SET);
        char *ret = fgets(l,
                          u_file.line_len*sizeof(char),
                          u_file.file);
        
        if (ret == NULL) {
            printf("NULL\n");
            return 1;
        }

        for (j = 0; j < u_file.num_fields; ++j)
            (*r)[j] = (*r)[j] | ((int)l[j*2] - 48);
    }

    return 0;
}
//}}}

//{{{int or_uint_records(struct uint_file u_file, 
int or_uint_fields(struct uint_file u_file, 
                    int num_f,
                    int *field_ids,
                    unsigned int **r)
{
    *r = (unsigned int *) calloc(u_file.num_records, sizeof(unsigned int));

    int a = 0;
    int i,j;
    char *pch, *ret;
    char *l = (char *) malloc(u_file.line_len*sizeof(char));

    for (i = 0; i < u_file.num_records; ++i) {
        //fprintf(stderr, "%d/%d\t",i,u_file.num_records);
        char *ret = fgets(l,
                          u_file.line_len*sizeof(char),
                          u_file.file);
        
        if (ret == NULL) {
            printf("NULL\n");
            return 1;
        }


        for (j = 0; j < num_f; ++j) {
            //fprintf(stderr, "%d ",((int)l[field_ids[j]*2] - 48)); 
            (*r)[i] = (*r)[i] | ((int)l[field_ids[j]*2] - 48);
        }

        //fprintf(stderr, "\n");
    }

    return 0;
}
//}}}

//{{{ int get_ubin_genotypes(struct ubin_file u_file,
int get_ubin_genotypes(struct ubin_file u_file,
                       int num_gts,
                       int *record_ids,
                       int *field_ids,
                       int *gts)
{
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    int i;
    int field_int_id, field_offset;
    unsigned int target_int;
    for (i = 0; i < num_gts; ++i) {
        field_int_id = field_ids[i]/16;
        field_offset = field_ids[i] - field_int_id*16;

        unsigned int target_int;
        // skip to the int that contains the target field
        fseek(u_file.file, 8 + // skip the record and field size field
                    record_ids[i]*num_bytes_per_record + // skip to the reccord
                    field_int_id*4,// skip to the propper int 
                    SEEK_SET);
        fread(&target_int,sizeof(unsigned int),1,u_file.file);

        // shift the int over to get the target field
        gts[i] = (target_int >> (30 - field_offset*2)) & 3;
    }

    return 0;
}
//}}}

//{{{ int or_ubin_fields(struct ubin_file u_file, 
int or_ubin_fields(struct ubin_file u_file, 
                   int num_f,
                   int *field_ids,
                   unsigned int **r)
{
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    // create enough space in r, and set all the values to zero
    *r = (unsigned int *) calloc(u_file.num_records, sizeof(unsigned int));

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    //unsigned int c[num_ints_per_record];
    
    int i,j;
    unsigned int target_int;
    for (i = 0; i < u_file.num_records; ++i) {

        fread(c,sizeof(unsigned int),num_ints_per_record,u_file.file);

        for (j = 0; j < num_f; ++j) {
            int field_int_id = field_ids[i]/16;
            int field_offset = field_ids[i] - field_int_id*16;

            target_int = c[field_int_id];
            (*r)[i] = ((*r)[i] | (target_int >> (30 - field_offset*2))) & 3;
        }

    }

    return 0;




    /*

    int i,j;
    *r = (unsigned int *) malloc(u_file.num_records*sizeof(unsigned int));
    int *record_ids = (int *) malloc(num_f*sizeof(int));
    int *gts = (int *) malloc(num_f*sizeof(int));

   for (j = 0; j < num_f; ++j) {
    record_ids[j] = i;
   }

    for (i = 0; i < u_file.num_records; ++i) {

        get_ubin_genotypes(u_file,
                           num_f,
                           record_ids,
                           field_ids,
                           gts);


        (*r)[i]=0;
        for (j = 0; j < num_f; ++j) {
            (*r)[i]= (*r)[i] | gts[j];
        }
    }

    free(record_ids);

    return 0;
    */
}
//}}}

//{{{ int or_ubin_records(struct ubin_file u_file, 
int or_ubin_records(struct ubin_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r)
{
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    // create enough space in r, and set all the values to zero
    *r = (unsigned int *) calloc(num_ints_per_record, sizeof(unsigned int));

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    //unsigned int c[num_ints_per_record];
    
    int i,j;
    for (i = 0; i < num_r; ++i) {

        // skip to the target record and read in the full record
        fseek(u_file.file, 8 + // skip the record and field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_ints_per_record,u_file.file);

        //and it
        for (j = 0; j < num_ints_per_record; ++j) {
            //printf("%u\t%u\t", (*r)[j], c[j]); 
            (*r)[j] = (*r)[j] | c[j];
            //printf("%u\n", (*r)[j]); 
        }

    }

    return 0;
}
//}}}

//{{{ void parse_cmd_line_int_csv(int *I,
void parse_cmd_line_int_csv(unsigned int *I,
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
