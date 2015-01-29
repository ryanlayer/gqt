/**
 * @file wah.c
 * @Author Ryan Layer (ryan.layer@gmail.com)
 * @date May, 2014
 * @brief Functions for converting and opperation on WAH-encoded files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <math.h>
#include <limits.h>
#include "genotq.h"

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
    int r = fread(&wf.num_fields,sizeof(uint32_t),1,wf.file);
    r = fread(&wf.num_records,sizeof(uint32_t),1,wf.file);

    wf.record_offsets = (uint64_t *) 
            malloc(sizeof (uint64_t)*wf.num_records);

    uint32_t i;
    for (i = 0; i < wf.num_records; ++i)
        r = fread(&(wf.record_offsets[i]),sizeof(uint64_t),1,wf.file);


    wf.header_offset = ftell(wf.file);

    return wf;
}
//}}}

// operate
//{{{ void wah_or(struct wah_run *x,
uint32_t  wah_or(struct wah_run *x,
                     struct wah_run *y,
                     uint32_t **O)
{
    struct wah_active_word a;
    uint32_t num_words;
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

    *O = (uint32_t *) malloc(Z_len * sizeof(uint32_t));

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

//{{{ uint32_t wah_and(struct wah_run *x,
uint32_t wah_and(struct wah_run *x,
                     struct wah_run *y,
                     uint32_t **O)
{
    struct wah_active_word a;
    uint32_t num_words;
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

    *O = (uint32_t *) malloc(Z_len * sizeof(uint32_t));

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

//{{{ struct wah_run init_wah_run(uint32_t *words,
struct wah_run init_wah_run(uint32_t *words,
                            uint32_t len){
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

// compress
//{{{ int ints_to_wah16(uint32_t *I,
uint32_t ints_to_wah16(uint32_t *I,
                           int I_len,
                           uint32_t used_bits,
                           uint16_t **W)
{
    uint32_t W_len;
    uint16_t *O;
    // split the intput up int to 31-bit groups
    uint32_t O_len = map_from_32_bits_to_15_bits(I, I_len, used_bits, &O);

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

//{{{ int ints_to_wah(uint32_t *I,
uint32_t ints_to_wah(uint32_t *I,
                         int I_len,
                         uint32_t used_bits,
                         uint32_t **W)
{
    uint32_t W_len;
    uint32_t *O;
    // split the intput up int to 31-bit groups
    uint32_t O_len = map_from_32_bits_to_31_bits(I, I_len, used_bits, &O);

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
    *W = (uint32_t *) malloc(W_len * sizeof(uint32_t));

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

//{{{unsinged int map_from_32_bits_to_15_bits(uint32_t *I,
uint32_t map_from_32_bits_to_15_bits(uint32_t *I,
                                         int I_len,
                                         uint32_t used_bits,
                                         uint16_t **O)
{
    uint32_t int_i, bit_i, group_i, in_group_i;
    uint32_t O_len =  (used_bits + 15 - 1)/ 15;
    //uint32_t O_len =  (I_len*32 + 31 - 1)/ 31;

    *O = (uint16_t *) calloc(O_len, sizeof(uint16_t));


    bit_i = 1;
    group_i = 0;
    in_group_i = 0;

    for (int_i = 0; int_i < I_len; ++int_i) {
        for ( ;bit_i<=32*(int_i+1); ++bit_i) {
            uint32_t bit = (I[int_i] >> (32 - (bit_i%32))) & 1;
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

//{{{unsinged int map_from_32_bits_to_31_bits(uint32_t *I,
/* 
 * Take a list of 32 bit number and gives the list of 31-bit groups
 * (represented by 32-bit with left padding) ints 
 * Returns the number of 31-bit groups in O
 */
uint32_t map_from_32_bits_to_31_bits(uint32_t *I,
                                         int I_len,
                                         uint32_t used_bits,
                                         uint32_t **O)
{
    uint32_t int_i, bit_i, group_i, in_group_i;
    uint32_t O_len =  (used_bits + 31 - 1)/ 31;
    //uint32_t O_len =  (I_len*32 + 31 - 1)/ 31;

    *O = (uint32_t *) calloc(O_len, sizeof(uint32_t));


    bit_i = 1;
    group_i = 0;
    in_group_i = 0;

    for (int_i = 0; int_i < I_len; ++int_i) {
        for ( ;bit_i<=32*(int_i+1); ++bit_i) {
            uint32_t bit = (I[int_i] >> (32 - (bit_i%32))) & 1;
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
            free(n);
            return 0;
        // The value on the tail is a fill of all zeros that is not full
        } else if ( ((*A_tail)->value.value >= 0x8000) &&
                    ((*A_tail)->value.value < 0xC000) ) {
            (*A_tail)->value.value += 1;
            free(n);
            return 0;
        } else { // the zeros cannot be added to the last active word
            (*A_tail)->next = n;
            *A_tail = n;
            return 1;
        }
    } else if (a.value == 0x7FFF) { // all ones
        if ( (*A_tail)->value.value == a.value ) {
            (*A_tail)->value.value = 0xC002;
            free(n);
            return 0;
        } else if ( (*A_tail)->value.value >= 0xC000) {
            (*A_tail)->value.value += 1;
            free(n);
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
            free(n);
            return 0;
        // The value on the tail is a fill of all zeros that is not full
        } else if ( ((*A_tail)->value.value >= 0x80000000) &&
                    ((*A_tail)->value.value < 0xC0000000) ) {
            (*A_tail)->value.value += 1;
            free(n);
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
            free(n);
            return 0;
        } else if ( (*A_tail)->value.value >= 0xC0000000) {
            (*A_tail)->value.value += 1;
            free(n);
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
                     uint32_t fill_size)
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

// inflate
//{{{ uint32_t wah_to_ints(uint32_t *W,
uint32_t wah_to_ints(uint32_t *W,
                         uint32_t W_len,
                         uint32_t **O)
{

    uint32_t wah_i;
    uint32_t num_bits = 0;

    for (wah_i = 0; wah_i < W_len; ++wah_i) {
        if (W[wah_i] >> 31 == 1) 
            num_bits += 31 * (W[wah_i] & 0x3fffffff); // zero out the fill bits
        else
            num_bits += 31;
    }

    uint32_t num_ints = (num_bits + 32 - 1) / 32;
    *O = (uint32_t *) malloc (num_ints * sizeof(uint32_t));


    uint32_t num_words,
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


//old
#if 0
//{{{ uint32_t get_wah_record(struct wah_file wf,
uint32_t get_wah_record(struct wah_file wf,
                        uint32_t wah_record,
                        uint32_t **wah)
{
    // get the size of the WAH-encoded bitmap
    uint32_t wah_size = 0, wah_offset = 0;
    if ( wah_record == 0) {
        wah_size = wf.record_offsets[wah_record];
        wah_offset = wf.header_offset;
    } else {
        wah_size = wf.record_offsets[wah_record] - 
                   wf.record_offsets[wah_record - 1];

        wah_offset = wf.header_offset +
                     sizeof(uint32_t) * 
                        (wf.record_offsets[wah_record] - wah_size);
    }

    *wah = (uint32_t *) malloc(sizeof(uint32_t)*wah_size);
    fseek(wf.file, wah_offset, SEEK_SET);
    int r = fread(*wah,sizeof(uint32_t),wah_size,wf.file);

    return wah_size;
}
//}}}

//{{{ uint32_t print_wah(struct wah_file wf,
uint32_t print_wah(struct wah_file wf,
                       uint32_t *record_ids,
                       uint32_t num_r,
                       uint32_t format)
{
    uint32_t i,j,k,wah_size,printed_bits,to_print = num_r;
    uint32_t *wah = NULL;

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
        uint32_t *ints = NULL;
        uint32_t ints_size = wah_to_ints(wah,wah_size,&ints);
        printed_bits = 0;

        for (j = 0; j < ints_size; ++j) {
            if (j !=0)
                printf(" ");
            for (k = 0; k < 16; ++k) {
                uint32_t val = (ints[j] >> (30 - 2*k)) & 3;
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

//{{{ uint32_t print_by_name_wah(char *wahbm_file_name,
uint32_t print_by_name_wah(char *wahbm_file_name,
                               uint32_t *record_ids,
                               uint32_t num_r,
                               uint32_t format)
{
    struct wah_file wf = init_wah_file(wahbm_file_name);
    return print_wah(wf, record_ids, num_r, format);
}
//}}} 
#endif
