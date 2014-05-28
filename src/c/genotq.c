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

//{{{ void plt_to_ubin(char *in_file_name, char *out_file_name);
int plt_by_name_to_ubin(char *in_file_name, char *out_file_name)
{

    struct plt_file pf = init_plt_file(in_file_name);
    int r = plt_to_ubin(pf, out_file_name);
    fclose(pf.file);
    return r;
}

int plt_to_ubin(struct plt_file pf, char *out_file_name)
{

    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);
    int i,j,count;
    size_t len;
    ssize_t read;
    int g[16];
    unsigned int *packed_ints;
    int num_packed_ints;

    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First value is the number of fields per record
    fwrite(&(pf.num_fields), sizeof(int), 1, o_file);

    // Second value is the number of records
    fwrite(&(pf.num_records), sizeof(int), 1, o_file);


    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        int r = plt_line_to_packed_ints(line,
                                        pf.num_fields,
                                        &packed_ints,
                                        &num_packed_ints);

        for (j = 0; j < num_packed_ints; ++j)
            fwrite(&(packed_ints[j]), sizeof(unsigned int), 1, o_file);

        free(packed_ints);
    }

    fclose(o_file);

    return 0;
}
//}}}

//{{{ int plt_line_to_packed_ints(char *line,
int plt_line_to_packed_ints(char *line,
                             int num_fields, 
                             unsigned int **packed_ints,
                             int *len)
{
    *len = 1 + ((num_fields - 1) / 16);
    *packed_ints = (unsigned *) malloc((*len)*sizeof(unsigned int));
    int i, two_bit_count, pack_int_count = 0;

    int g[16];

    memset(g,0,sizeof(g));
    two_bit_count = 0;
    for (i = 0; i < num_fields; ++i) {
        g[i%16] = ((int)line[i*2] - 48);

        two_bit_count +=1;

        if (two_bit_count == 16) {
            (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);
            pack_int_count += 1;
            memset(g,0,sizeof(g));
            two_bit_count = 0;
        }
        
    }

    if (two_bit_count > 0)
        (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);

    return 0;
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
        fseek(uf.file, uf.header_offset + // skip the record and field size field
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

//{{{unsinged int map_from_32_bits_to_31_bits(unsigned int *I,
/* 
 * Take a list of 32 bit number and gives the list of 31-bit groups
 * (represented by 32-bit with left padding) ints 
 * Returns the number of 31-bit groups in O
 */
unsigned int map_from_32_bits_to_31_bits(unsigned int *I,
                                int I_len,
                                unsigned int **O)
{
    unsigned int int_i, bit_i, group_i, in_group_i;
    unsigned int total_bits = I_len *32;
    unsigned int O_len =  (total_bits + 31 - 1)/ 31;

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
        }
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

//{{{ int ints_to_wah(unsigned int *I,
unsigned int ints_to_wah(unsigned int *I,
                 int I_len,
                 unsigned int **W)
{
    unsigned int W_len;
    unsigned int *O;
    // split the intput up int to 31-bit groups
    unsigned int O_len = map_from_32_bits_to_31_bits(I, I_len, &O);

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
            wah_run_decode(x);
        }

        // IF (yrun.nWords == 0) ++yrun.it, yrun.decode();
        if (y->num_words == 0) {
            y->word_i += 1;
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
void parse_cmd_line_int_csv(int *I,
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