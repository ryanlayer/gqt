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
#include <sysexits.h>
#include "genotq.h"
#ifdef __SSE4_2__
#include <nmmintrin.h>
#endif

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
    uint32_t rle_v = 0;
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
                if (!n)
                    err(EX_OSERR, "malloc error");
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
    if (!n)
        err(EX_OSERR, "malloc error");
    n->value = rle_v;
    n->next = NULL;
    ll_len += 1;

    if (head == NULL)
        head = n;
    else
        tail->next = n;

    tail = n;

    *O = (uint32_t *) malloc(ll_len*sizeof(uint32_t));
    if (!O)
        err(EX_OSERR, "malloc error");
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
    if (!r)
        err(EX_OSERR, "malloc error");

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
    if (!r)
        err(EX_OSERR, "malloc error");

    int i;
    for (i = 0; i < 16; ++i) 
        r[i] = (packed_ints >> (30 - i*2)) & 3;
    
    return r;
}
//}}}

//{{{int popcount(uint32_t x) {
int popcount(uint32_t x) {

#ifndef __SSE4_2__
    return __builtin_popcount(x);
    /*
    int count;
    for (count=0; x; count++)
        x &= x-1;
    return count;
    */
#else
    return _mm_popcnt_u32(x);
#endif

}
//}}}

//{{{ void check_file_read(char *file_name, FILE *fp, size_t exp, size_t obs)
void check_file_read(char *file_name, FILE *fp, size_t exp, size_t obs)
{
    if (exp != obs) {
        if (feof(fp))
            errx(EX_IOERR,
                 "Error reading file \"%s\": End of file",
                 file_name);
        err(EX_IOERR, "Error reading file \"%s\"", file_name);
    }
}
//}}}

//{{{ int check_field_name(char *field_name)
int check_field_name(char *field_name)
{
    // The first character cannot be a numer

    if ((field_name[0] >= '0') && (field_name[0] <= '9'))
        return 0;

    int i;

    for (i = 0; i < strlen(field_name); ++i) {
        if ( (field_name[i] < '0') ||
            ((field_name[i] >= ':') && (field_name[i] <= '@')) ||
            ((field_name[i] >= '[') && (field_name[i] <= '`') &&
                (field_name[i] != '_')) ||
             (field_name[i] > 'z') )
            return i;
    }

    return -1;
}
//}}}

//{{{ int is_int(char *s, int *v)
//base on http://rus.har.mn/blog/2014-05-19/strtol-error-checking/
// 1: is an int
// 0: is text
int is_int(char *s, int *v)
{
    errno = 0;
    char *endptr;
    long val = strtol(s, &endptr, 10);
    if ( ((errno != 0 ) ||(*endptr != '\0')) || (val>INT_MAX))
        return 0;
    else {
        *v = (int) val;
        return 1;
    }
}
//}}}

//{{{struct gqt_file_header new_gqt_file_header(char type)
struct gqt_file_header *new_gqt_file_header(char type,
                                            char *full_cmd,
                                            uint32_t num_variants,
                                            uint32_t num_samples)
{
    struct gqt_file_header *h = (struct gqt_file_header *) 
            malloc(sizeof(struct gqt_file_header));

    h->marker[0] = 'G';
    h->marker[1] = 'Q';
    h->marker[2] = 'T';

    if ( !((type != 'g') || (type != 'v') || (type != 'b')) )
        errx(EX_UNAVAILABLE,
             "Cannot create header for unknown file type '%c'.",
             type);

    h->type = type;

    h->major = atoi(MAJOR_VERSION);
    h->minor = atoi(MINOR_VERSION);
    h->revision = atoi(REVISION_VERSION);
    h->build = atoi(BUILD_VERSION);
    h->magic  = 0x11223344;
    h->num_variants = num_variants;
    h->num_samples = num_samples;
    h->id_hash = hash_cmd(full_cmd);
    memset(h->more, 0, MORE_SIZE*sizeof(uint32_t));

    return h;
}
//}}}

//{{{struct gqt_file_header read_gqt_file_header(FILE *f)
struct gqt_file_header *read_gqt_file_header(char *file_name, FILE *f)
{
    struct gqt_file_header *h = (struct gqt_file_header *) 
            malloc(sizeof(struct gqt_file_header));

    if (fseek(f, 0, SEEK_SET))
        err(EX_IOERR, "Error seeking to header in VID file '%s'.", file_name);

    size_t fr = fread(h, sizeof(struct gqt_file_header), 1, f);
    check_file_read(file_name, f, 1, fr);
    return h;
}
//}}}

//{{{ char *hash_cmd(char *full_cmd)
unsigned long hash_cmd(char *full_cmd)
{
    //djb2
    unsigned long hash = 5381;
    int c;

    while ((c = *full_cmd++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}
//}}}
