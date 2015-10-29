/**
 * @file genotq.c
 * @Author Ryan Layer (ryan.layer@gmail.com)
 * @date May, 2014
 * @brief Functions for converting and opperation on various encoding
 * strategies
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>
#include <assert.h>
#include <immintrin.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>

#include "wahbm.h"
#include "wahbm_in_place.h"
#include "genotq.h"
#include "pthread_pool.h"
#include "timer.h"

const int tab32[32] = {
    0,  9,  1, 10, 13, 21,  2, 29,
    11, 14, 16, 18, 22, 25,  3, 30,
    8, 12, 20, 28, 15, 17, 24,  7,
    19, 27, 23,  6, 26,  5,  4, 31};

//{{{ int log2_32 (uint32_t value)
int log2_32 (uint32_t value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return tab32[(uint32_t)(value*0x07C4ACDD) >> 27];
}
//}}}

// wahbm
//{{{ struct wahbm_file *new_wahbm_file(char *file_name,
struct wahbm_file *new_wahbm_file(char *file_name,
                                  char *full_cmd,
                                  uint32_t num_variants,
                                  uint32_t num_samples)
{
    struct wahbm_file *wf = (struct wahbm_file *) 
            malloc(sizeof(struct wahbm_file));

    wf->type = WAHBM_LOCAL;
    wf->file_name = strdup(file_name);

    wf->file.local = fopen(file_name,"wb");
    if (!wf->file.local)
        err(EX_CANTCREAT, "Cannot create WAHBM file '%s'", file_name);

    wf->gqt_header = new_gqt_file_header('w', 
                                         full_cmd,
                                         num_variants,
                                         num_samples);

    //wf->word_size = 31;

    if (fwrite(wf->gqt_header,
               sizeof(struct gqt_file_header),
               1,
               wf->file.local) != 1)
        err(EX_IOERR, "Error writing header to WAHBM file '%s'", file_name);

    // 4 BM for each individual
    wf->record_offsets = 
            (uint64_t *) calloc(4*num_samples, sizeof(uint64_t));


    if (fwrite(wf->record_offsets,
               sizeof(uint64_t),
               4*num_samples,
               wf->file.local) != 4*num_samples)
        err(EX_IOERR, "Error writing header to WAHBM file '%s'", file_name);

    wf->header_offset = ftell(wf->file.local);

    return wf;
}
//}}}

//{{{void append_ubin_to_wahbm_file(struct wahbm_file *wf,
void append_ubin_to_wahbm_file(struct wahbm_file *wf,
                               uint32_t wah_i,
                               uint32_t *ubin,
                               uint32_t ubin_len)
{
    assert(wf->type == WAHBM_LOCAL);

    uint32_t *wah;
    uint32_t *wah_sizes;
    uint32_t wah_len = ubin_to_bitmap_wah(ubin,
                                          ubin_len,
                                          wf->gqt_header->num_variants,
                                          &wah,
                                          &wah_sizes);


    uint64_t offset_total = 0;
    if (wah_i > 0)
        offset_total = wf->record_offsets[wah_i*4 - 1];

    uint32_t i;
    
    // Update index 
    for (i = 0; i < 4; ++i) {
        offset_total += wah_sizes[i];
        wf->record_offsets[wah_i*4 + i] = offset_total;
    }

    // Jump to the end of the file
    if( fseek(wf->file.local, 0, SEEK_END) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            wf->file_name);

    // Write out the compressed WAH bitmaps
    if (fwrite(wah, 
               sizeof(uint32_t), 
               wah_len,
               wf->file.local) != wah_len)
        err(EX_IOERR, "Error writing WAH to '%s'.", wf->file_name); 

    free(wah);
    free(wah_sizes);
}
//}}}

//{{{void append_wahbm_to_wahbm_file(struct wahbm_file *wf,
void append_wahbm_to_wahbm_file(struct wahbm_file *wf,
                                uint32_t wah_i,
                                uint32_t *wah,
                                uint32_t wah_len,
                                uint32_t *wah_sizes)
{
    assert(wf->type == WAHBM_LOCAL);

    uint64_t offset_total = 0;
    if (wah_i > 0)
        offset_total = wf->record_offsets[wah_i*4 - 1];

    uint32_t i;
    
    // Update index 
    for (i = 0; i < 4; ++i) {
        offset_total += wah_sizes[i];
        wf->record_offsets[wah_i*4 + i] = offset_total;
    }

    // Jump to the end of the file
    if( fseek(wf->file.local, 0, SEEK_END) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            wf->file_name);

    // Write out the compressed WAH bitmaps
    if (fwrite(wah, 
               sizeof(uint32_t), 
               wah_len,
               wf->file.local) != wah_len)
        err(EX_IOERR, "Error writing WAH to '%s'.", wf->file_name); 
}
//}}}

//{{{void write_record_offsets(struct wahbm_file *wf)
void write_record_offsets(struct wahbm_file *wf)
{
    assert(wf->type == WAHBM_LOCAL);

    // Jump to the correct point in the WAH header
    if ( fseek(wf->file.local, sizeof(struct gqt_file_header), SEEK_SET) )
        err(EX_IOERR,
            "Error seeking to record offsets in WAH BM file '%s'.",
            wf->file_name);

    if (fwrite(wf->record_offsets,
               sizeof(uint64_t),
               4*(wf->gqt_header->num_samples),
               wf->file.local) != 4*(wf->gqt_header->num_samples))
        err(EX_IOERR,
            "Error writing header to WAHBM file '%s'",
            wf->file_name);
}
//}}}

//{{{void destroy_wahbm_file(struct wahbm_file *wf)
void destroy_wahbm_file(struct wahbm_file *wf)
{
    if (wf->type == WAHBM_LOCAL) {
        if (fclose(wf->file.local))
            err(EX_IOERR, "Error closeing WAHBM file '%s'", wf->file_name); 
    } else {
        if (knet_close(wf->file.remote) != 0)
            err(EX_IOERR,
                "Error closeing remote WAHBM file '%s'",
                wf->file_name); 
    }
    free(wf->record_offsets);
    free(wf);
}
//}}}

//{{{ struct wahbm_file *open_wahbm_file(char *file_name)
struct wahbm_file *open_wahbm_file(char *file_name)
{
    struct wahbm_file *wf = (struct wahbm_file *) 
            malloc(sizeof(struct wahbm_file));

    wf->type = WAHBM_REMOTE;
    wf->file_name = strdup(file_name);
    wf->file.remote = knet_open(file_name,"rb+");

    if (!(wf->file.remote))
        err(EX_NOINPUT, "Cannot open WAHBM file '%s'", file_name);

    wf->gqt_header = read_remote_gqt_file_header(wf->file_name,
                                                 wf->file.remote);

    if ( !((wf->gqt_header->marker[0] == 'G') &&
           (wf->gqt_header->marker[1] == 'Q') && 
           (wf->gqt_header->marker[2] == 'T')) )
        errx(EX_NOINPUT, "File '%s' is not a GQT file.", file_name);

    if (wf->gqt_header->type != 'w')
        errx(EX_NOINPUT, "File '%s' is not a WAHBM file.", file_name);

    // 4 BM for each individual
    wf->record_offsets = (uint64_t *) 
            malloc(4*wf->gqt_header->num_samples*sizeof(uint64_t));

    size_t fr = knet_read(wf->file.remote,
                          wf->record_offsets,
                          4*wf->gqt_header->num_samples * sizeof(uint64_t));

    check_remote_file_read(wf->file_name,
                           4*wf->gqt_header->num_samples * sizeof(uint64_t),
                           fr);

    /*
    size_t fr = fread(wf->record_offsets,
                      sizeof(uint64_t),
                      4*wf->gqt_header->num_samples,
                      wf->file);
    check_file_read(wf->file_name, wf->file, 4*wf->gqt_header->num_samples, fr);
    */

    //wf->header_offset = ftell(wf->file);
    wf->header_offset = knet_tell(wf->file.remote);

    return wf;
}
//}}}

//{{{ uint32_t get_wahbm_bitmap(struct wfile *wf,
uint32_t get_wahbm_bitmap(struct wahbm_file *wf,
                          uint32_t bitmap_i,
                          uint32_t bitarray_i,
                          uint32_t **wah_bitmap)
{

    if (*wah_bitmap == NULL) {
        uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
        *wah_bitmap = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
        if (!(*wah_bitmap))
            err(EX_OSERR, "malloc error");
    }

    // get the size of the WAH-encoded bitmap
    uint64_t wah_size = 0;
    uint64_t wah_offset = 0;
    if ((bitmap_i == 0) && (bitarray_i == 0)) {
        wah_size = wf->record_offsets[bitmap_i + bitarray_i];
        wah_offset = wf->header_offset;
    } else {
        wah_size = wf->record_offsets[bitmap_i*4 + bitarray_i] - 
                wf->record_offsets[bitmap_i*4 + bitarray_i - 1];

        wah_offset = wf->header_offset + sizeof(uint32_t) * 
                (wf->record_offsets[bitmap_i*4 + bitarray_i] - wah_size);
    }

    if (wf->type == WAHBM_LOCAL) {
        fseek(wf->file.local, wah_offset, SEEK_SET);
        size_t fr = fread(*wah_bitmap,sizeof(uint32_t),wah_size,wf->file.local);
        check_file_read(wf->file_name, wf->file.local, wah_size, fr);
    } else {
        if (knet_seek(wf->file.remote, wah_offset, SEEK_SET) == -1)
            err(EX_IOERR,
                "Error seeking to data in remote WAHBM file '%s'.",
                wf->file_name);

        size_t fr = knet_read(wf->file.remote,
                              *wah_bitmap,
                              wah_size*sizeof(uint32_t));
        check_remote_file_read(wf->file_name,
                               wah_size*sizeof(uint32_t),
                               fr);
    }

    return (uint32_t)wah_size;
}
//}}}

#if 0
//{{{ struct wah_file init_wahbm_file(char *file_name)
struct wah_file init_wahbm_file(char *file_name)
{
    struct wah_file wf;

    wf.file_name = strdup(file_name);

    wf.file = fopen(file_name, "rb");
    if (!wf.file)
        err(EX_NOINPUT, "Cannot open file \"%s\"", file_name);


    // Jump to the begining of the file to grab the record size
    fseek(wf.file, 0, SEEK_SET);
    size_t fr = fread(&wf.num_fields,sizeof(uint32_t),1,wf.file);
    check_file_read(file_name, wf.file, 1, fr);

    fr = fread(wf.num_records,sizeof(uint32_t),1,wf.file);
    check_file_read(file_name, wf.file, 1, fr);

    wf.record_offsets = (uint64_t *) 
            malloc(sizeof (uint64_t)*wf.num_records*4);
    if (!wf.record_offsets)
        err(EX_OSERR, "malloc error");

    uint32_t i;
    for (i = 0; i < wf.num_records*4; ++i) {
        fr = fread(&(wf.record_offsets[i]),sizeof(uint64_t),1,wf.file);
        check_file_read(file_name, wf.file, 1, fr);
    }


    wf.header_offset = ftell(wf.file);

    return wf;
}
//}}}
#endif

//{{{ uint32_t wahbm_speed_check(char *in)
uint32_t wahbm_speed_check(char *in)
{
    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_0_s, i_1_s, i_2_s;
    uint32_t *i_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_0)
        err(EX_OSERR, "malloc error");
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1)
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2)
        err(EX_OSERR, "malloc error");

    uint32_t j_0_s, j_1_s, j_2_s;
    uint32_t *j_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_0)
        err(EX_OSERR, "malloc error");
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1)
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2)
        err(EX_OSERR, "malloc error");

    uint32_t x_0_s, x_1_s, x_2_s;
    uint32_t *x_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_0)
        err(EX_OSERR, "malloc error");
    uint32_t *x_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1)
        err(EX_OSERR, "malloc error");
    uint32_t *x_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2)
        err(EX_OSERR, "malloc error");


    uint32_t x_0_sc, x_1_sc, x_2_sc;
    uint32_t *x_0_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_0_c)
        err(EX_OSERR, "malloc error");
    uint32_t *x_1_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1_c)
        err(EX_OSERR, "malloc error");
    uint32_t *x_2_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2_c)
        err(EX_OSERR, "malloc error");

    uint32_t t_s;
    uint32_t *t = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!t)
        err(EX_OSERR, "malloc error");

    uint32_t i,j,k, i_to_j_d;

        
    //i_0_s = get_wah_bitmap_in_place(wf, 0, 0, &i_0);
    i_0_s = get_wahbm_bitmap(wf, 0, 0, &i_0);
    //i_1_s = get_wah_bitmap_in_place(wf, 0, 1, &i_1);
    i_1_s = get_wahbm_bitmap(wf, 0, 1, &i_1);
    //i_2_s = get_wah_bitmap_in_place(wf, 0, 2, &i_2);
    i_2_s = get_wahbm_bitmap(wf, 0, 2, &i_2);

    memset(x_0_c, 0, sizeof(uint32_t)*max_wah_size);
    memset(x_1_c, 0, sizeof(uint32_t)*max_wah_size);
    memset(x_2_c, 0, sizeof(uint32_t)*max_wah_size);


    x_0_sc = wah_in_place_or(x_0_c, max_wah_size, i_0, i_0_s);
    x_1_sc = wah_in_place_or(x_1_c, max_wah_size, i_1, i_1_s);
    x_2_sc = wah_in_place_or(x_2_c, max_wah_size, i_2, i_2_s);

    printf("\nXOR\n");
    for (i = 0; i < 5; ++i) {
    start();
    for (j = 1; j < wf->gqt_header->num_variants; ++j) {

        memcpy(x_0, x_0_c, x_0_sc * sizeof(uint32_t));
        memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
        memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
        x_0_s = x_0_sc;
        x_1_s = x_1_sc;
        x_2_s = x_2_sc;

        //load in the jth record
        //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
        j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
        //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
        j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
        //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
        j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
        //find the xor of all 4
        
        x_0_s =  wah_in_place_xor(x_0, x_0_s, j_0, j_0_s);
        x_1_s =  wah_in_place_xor(x_1, x_1_s, j_1, j_1_s);
        x_2_s =  wah_in_place_xor(x_2, x_2_s, j_2, j_2_s);
    }
    stop();
    printf("%lu\n", report());
    }

    printf("\nAND\n");
    for (i = 0; i < 5; ++i) {
    start();
    for (j = 1; j < wf->gqt_header->num_variants; ++j) {

        memcpy(x_0, x_0_c, x_0_sc * sizeof(uint32_t));
        memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
        memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
        x_0_s = x_0_sc;
        x_1_s = x_1_sc;
        x_2_s = x_2_sc;

        //load in the jth record
        //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
        j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
        //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
        j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
        //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
        j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
        //find the xor of all 4
        
        x_0_s =  wah_in_place_and(x_0, x_0_s, j_0, j_0_s);
        x_1_s =  wah_in_place_and(x_1, x_1_s, j_1, j_1_s);
        x_2_s =  wah_in_place_and(x_2, x_2_s, j_2, j_2_s);
    }
    stop();
    printf("%lu\n", report());
    }

    printf("\nOR\n");
    for (i = 0; i < 5; ++i) {
    start();
    for (j = 1; j < wf->gqt_header->num_variants; ++j) {

        memcpy(x_0, x_0_c, x_0_sc * sizeof(uint32_t));
        memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
        memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
        x_0_s = x_0_sc;
        x_1_s = x_1_sc;
        x_2_s = x_2_sc;

        //load in the jth record
        //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
        j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
        //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
        j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
        //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
        j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
        //find the xor of all 4
        
        x_0_s =  wah_in_place_or(x_0, x_0_s, j_0, j_0_s);
        x_1_s =  wah_in_place_or(x_1, x_1_s, j_1, j_1_s);
        x_2_s =  wah_in_place_or(x_2, x_2_s, j_2, j_2_s);
    }
    stop();
    printf("%lu\n", report());
    }

    printf("\nOR\n");
    for (i = 0; i < 5; ++i) {
    start();
    for (j = 1; j < wf->gqt_header->num_variants; ++j) {

        //memcpy(x_0, x_0_c, x_0_sc * sizeof(uint32_t));
        //memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
        //memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
        //x_0_s = x_0_sc;
        //x_1_s = x_1_sc;
        //x_2_s = x_2_sc;

        //load in the jth record
        //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
        j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
        //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
        j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
        //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
        j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
        //find the xor of all 4
        
        x_0_s =  wah_or_b(x_0, i_0, i_0_s, j_0, j_0_s);
        x_1_s =  wah_or_b(x_1, i_1, i_1_s, j_1, j_1_s);
        x_2_s =  wah_or_b(x_2, i_2, i_2_s, j_2, j_2_s);
    }
    stop();
    printf("%lu\n", report());
    }

    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}
//}}}

//{{{ uint32_t wahbm_pca_by_name(char *in, char *out)
uint32_t wahbm_pca_by_name(char *in, char *out)
{

    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);

    //uint32_t max_wah_size = (wf.num_fields + 31 - 1)/ 31;
    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_0_s, i_1_s, i_2_s;
    uint32_t *i_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_0 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2 )
        err(EX_OSERR, "malloc error");

    uint32_t j_0_s, j_1_s, j_2_s;
    uint32_t *j_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_0 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2 )
        err(EX_OSERR, "malloc error");

    uint32_t x_0_s, x_1_s, x_2_s;
    uint32_t *x_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_0 )
        err(EX_OSERR, "malloc error");
    uint32_t *x_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *x_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2 )
        err(EX_OSERR, "malloc error");


    uint32_t x_0_sc, x_1_sc, x_2_sc;
    uint32_t *x_0_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_0_c )
        err(EX_OSERR, "malloc error");
    uint32_t *x_1_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1_c )
        err(EX_OSERR, "malloc error");
    uint32_t *x_2_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2_c )
        err(EX_OSERR, "malloc error");

    uint32_t t_s;
    uint32_t *t = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!t)
        err(EX_OSERR, "malloc error");

    uint32_t i,j,k, i_to_j_d;

#ifdef time_wahbm_pca_by_name
    unsigned long t_i_get_wah_bitmap_in_place   = 0,
                  t_i_memset                    = 0,
                  t_i_wah_in_place_or           = 0,
                  t_j_memcpy                    = 0,
                  t_j_get_wah_bitmap_in_place   = 0,
                  t_j_wah_in_place_xor          = 0,
                  t_j_r1                        = 0,
                  t_j_r2                        = 0,
                  t_j_popcounts                 = 0,
                  t_j_print                     = 0;
#endif

    //for (i = 0; i < wf.num_records; ++i) {
    for (i = 0; i < wf->gqt_header->num_variants; ++i) {
        //load in the ith record

#ifdef time_wahbm_pca_by_name
        start();
#endif

        //i_0_s = get_wah_bitmap_in_place(wf, i, 0, &i_0);
        i_0_s = get_wahbm_bitmap(wf, i, 0, &i_0);
        //i_1_s = get_wah_bitmap_in_place(wf, i, 1, &i_1);
        i_1_s = get_wahbm_bitmap(wf, i, 1, &i_1);
        //i_2_s = get_wah_bitmap_in_place(wf, i, 2, &i_2);
        i_2_s = get_wahbm_bitmap(wf, i, 2, &i_2);

#ifdef time_wahbm_pca_by_name
        stop();
        t_i_get_wah_bitmap_in_place += report();
#endif


#ifdef time_wahbm_pca_by_name
        start();
#endif

        memset(x_0_c, 0, sizeof(uint32_t)*max_wah_size);
        memset(x_1_c, 0, sizeof(uint32_t)*max_wah_size);
        memset(x_2_c, 0, sizeof(uint32_t)*max_wah_size);

#ifdef time_wahbm_pca_by_name
        stop();
        t_i_memset += report();
#endif

#ifdef time_wahbm_pca_by_name
        start();
#endif

        x_0_sc = wah_in_place_or(x_0_c, max_wah_size, i_0, i_0_s);
        x_1_sc = wah_in_place_or(x_1_c, max_wah_size, i_1, i_1_s);
        x_2_sc = wah_in_place_or(x_2_c, max_wah_size, i_2, i_2_s);

#ifdef time_wahbm_pca_by_name
        stop();
        t_i_wah_in_place_or += report();
#endif

        //for (j = i+1; j < wf.num_records; ++j) {
        for (j = i+1; j < wf->gqt_header->num_variants; ++j) {

#ifdef time_wahbm_pca_by_name
            start();
#endif

            memcpy(x_0, x_0_c, x_0_sc * sizeof(uint32_t));
            memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
            memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
            x_0_s = x_0_sc;
            x_1_s = x_1_sc;
            x_2_s = x_2_sc;

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_memcpy += report();
#endif

#ifdef time_wahbm_pca_by_name
            start();
#endif

            //load in the jth record
            //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
            j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
            //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
            j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
            //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
            j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
            //find the xor of all 4
            
#ifdef time_wahbm_pca_by_name
            stop();
            t_j_get_wah_bitmap_in_place += report();
#endif

#ifdef time_wahbm_pca_by_name
            start();
#endif

            x_0_s =  wah_in_place_xor(x_0, x_0_s, j_0, j_0_s);
            x_1_s =  wah_in_place_xor(x_1, x_1_s, j_1, j_1_s);
            x_2_s =  wah_in_place_xor(x_2, x_2_s, j_2, j_2_s);

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_wah_in_place_xor += report();
#endif

            memcpy(t, x_0, x_0_s * sizeof(uint32_t));
            t_s = x_0_s;

#ifdef time_wahbm_pca_by_name
            start();
#endif

            //r1-> (0 AND 1) OR (1 AND 2) --> (0 OR 2) AND 1
            x_0_s =  wah_in_place_or(x_0, x_0_s, x_2, x_2_s);
            x_0_s =  wah_in_place_and(x_0, x_0_s, x_1, x_1_s);

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_r1 += report();
#endif
            
#ifdef time_wahbm_pca_by_name
            start();
#endif
            //r2-> 2 AND 0
            x_2_s =  wah_in_place_and(x_2, x_2_s, t, t_s);

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_r2 += report();
#endif

#ifdef time_wahbm_pca_by_name
            start();
#endif

            i_to_j_d = 0;
            /*
            for (k = 0; k < x_0_s; ++k)
                i_to_j_d += popcount(x_0[k]);
            for (k = 0; k < x_2_s; ++k)
                i_to_j_d += popcount(x_2[k])*4;
            */

            uint64_t *l_x_0 = (uint64_t *)x_0;
            uint64_t *l_x_2 = (uint64_t *)x_2;

            for (k = 0; k < x_0_s/2; k++)
                i_to_j_d += __builtin_popcountll(l_x_0[k]);
            if (k*2 < x_0_s)
                i_to_j_d += __builtin_popcount(x_0[k*2]);

            for (k = 0; k < x_2_s/2; k++)
                i_to_j_d += __builtin_popcountll(l_x_2[k])*4;
            if (k*2 < x_2_s)
                i_to_j_d += __builtin_popcount(x_2[k*2])*4;

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_popcounts += report();
#endif

#ifdef time_wahbm_pca_by_name
            start();
#endif

            if (j!=i+1)
                printf("\t");
            //printf("%f", ((float)i_to_j_d)/((float)wf.num_fields));
            printf("%f",
                   ((float)i_to_j_d)/((float)wf->gqt_header->num_variants));

#ifdef time_wahbm_pca_by_name
            stop();
            t_j_print += report();
#endif
        }
        printf("\n");
    }

#ifdef time_wahbm_pca_by_name
    float t_total = t_i_get_wah_bitmap_in_place   +
                    t_i_memset                    +
                    t_i_wah_in_place_or           +
                    t_j_memcpy                    +
                    t_j_get_wah_bitmap_in_place   +
                    t_j_wah_in_place_xor          +
                    t_j_r1                        +
                    t_j_r2                        +
                    t_j_popcounts                 +
                    t_j_print;

    fprintf(stderr, "t_i_get_wah_bitmap_in_place:%lu\t%f\n"
                    "t_i_memset:%lu\t%f\n"
                    "t_i_wah_in_place_or:%lu\t%f\n"
                    "t_j_memcpy:%lu\t%f\n"
                    "t_j_get_wah_bitmap_in:%lu\t%f\n"
                    "t_j_wah_in_place_xor:%lu\t%f\n"
                    "t_j_r1:%lu\t%f\n"
                    "t_j_r2:%lu\t%f\n"
                    "t_j_popcounts:%lu\t%f\n"
                    "t_j_print:%lu\t%f\n",
                    t_i_get_wah_bitmap_in_place,
                    t_i_get_wah_bitmap_in_place/t_total,
                    t_i_memset,
                    t_i_memset/t_total,
                    t_i_wah_in_place_or,
                    t_i_wah_in_place_or/t_total,
                    t_j_memcpy,
                    t_j_memcpy/t_total,
                    t_j_get_wah_bitmap_in_place,
                    t_j_get_wah_bitmap_in_place/t_total,
                    t_j_wah_in_place_xor,
                    t_j_wah_in_place_xor/t_total,
                    t_j_r1,
                    t_j_r1/t_total,
                    t_j_r2,
                    t_j_r2/t_total,
                    t_j_popcounts,
                    t_j_popcounts/t_total,
                    t_j_print,
                    t_j_print/t_total);
#endif


    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}
//}}}

//{{{ uint32_t wahbm_hamm_dist_by_name(char *in, char *out)
uint32_t wahbm_hamm_dist_by_name(char *in, char *out)
{

    FILE *f = fopen(out, "w");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", out);

    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);

    //uint32_t max_wah_size = (wf.num_fields + 31 - 1)/ 31;
    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_0_s, i_1_s, i_2_s, i_a_s;
    uint32_t *i_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_0 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2 )
        err(EX_OSERR, "malloc error");

    uint32_t *i_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_a )
        err(EX_OSERR, "malloc error");

    uint32_t j_0_s, j_1_s, j_2_s;
    uint32_t *j_0 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_0 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2 )
        err(EX_OSERR, "malloc error");

    uint32_t j_r_s, j_a_s, j_het_s, j_hom_s;
    uint32_t *j_r = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_r )
        err(EX_OSERR, "malloc error");
    uint32_t *j_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_a )
        err(EX_OSERR, "malloc error");
    uint32_t *j_het = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_het )
        err(EX_OSERR, "malloc error");
    uint32_t *j_hom = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_hom )
        err(EX_OSERR, "malloc error");


    uint32_t i,j,k, i_to_j_d;

    //for (i = 0; i < wf.num_records; ++i) {
    for (i = 0; i < wf->gqt_header->num_variants; ++i) {
        //load in the ith record

        //i_0_s = get_wah_bitmap_in_place(wf, i, 0, &i_0);
        i_0_s = get_wahbm_bitmap(wf, i, 0, &i_0);
        //i_1_s = get_wah_bitmap_in_place(wf, i, 1, &i_1);
        i_1_s = get_wahbm_bitmap(wf, i, 1, &i_1);
        //i_2_s = get_wah_bitmap_in_place(wf, i, 2, &i_2);
        i_2_s = get_wahbm_bitmap(wf, i, 2, &i_2);

        memset(i_a, 0, sizeof(uint32_t)*max_wah_size);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_1, i_1_s);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_2, i_2_s);

        //for (j = i+1; j < wf.num_records; ++j) {
        for (j = i+1; j < wf->gqt_header->num_variants; ++j) {
            //load in the jth record
            //j_0_s = get_wah_bitmap_in_place(wf, j, 0, &j_0);
            j_0_s = get_wahbm_bitmap(wf, j, 0, &j_0);
            //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
            j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
            //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
            j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);


            memset(j_r,   0, sizeof(uint32_t)*max_wah_size);
            memset(j_a,   0, sizeof(uint32_t)*max_wah_size);
            memset(j_het, 0, sizeof(uint32_t)*max_wah_size);
            memset(j_hom, 0, sizeof(uint32_t)*max_wah_size);

            j_r_s = wah_in_place_or(j_r, max_wah_size, j_0, j_0_s);

            j_a_s = wah_in_place_or(j_a, max_wah_size, j_1, j_1_s);
            j_a_s = wah_in_place_or(j_a, max_wah_size, j_2, j_2_s);

            j_het_s = wah_in_place_or(j_het, max_wah_size, j_1, j_1_s);

            j_hom_s = wah_in_place_or(j_hom, max_wah_size, j_2, j_2_s);

            // 0 -> 1/2
            j_a_s = wah_in_place_and(j_a, max_wah_size, i_0, i_0_s);

            // 1/2 -> 0
            j_r_s = wah_in_place_and(j_r, max_wah_size, i_a, i_a_s);

            // 2 -> 1
            j_het_s = wah_in_place_and(j_het, max_wah_size, i_2, i_2_s);
            
            // 1 -> 2
            j_hom_s = wah_in_place_and(j_hom, max_wah_size, i_1, i_1_s);

            i_to_j_d = 0;

            /*
            for (k = 0; k < max_wah_size; ++k) {
                i_to_j_d += popcount(j_a[k]);
                i_to_j_d += popcount(j_r[k]);
                i_to_j_d += popcount(j_het[k]);
                i_to_j_d += popcount(j_hom[k]);
            }
            */

            uint64_t *l_j_a = (uint64_t *)j_a;
            uint64_t *l_j_r = (uint64_t *)j_r;
            uint64_t *l_j_het = (uint64_t *)j_het;
            uint64_t *l_j_hom = (uint64_t *)j_hom;

            for (k = 0; k < j_a_s/2; k++) {
                i_to_j_d += __builtin_popcountll(l_j_a[k]);
                i_to_j_d += __builtin_popcountll(l_j_r[k]);
                i_to_j_d += __builtin_popcountll(l_j_het[k]);
                i_to_j_d += __builtin_popcountll(l_j_hom[k]);
            }

            if (k*2 < j_r_s) {
                i_to_j_d += __builtin_popcount(j_a[k*2]);
                i_to_j_d += __builtin_popcount(j_r[k*2]);
                i_to_j_d += __builtin_popcount(j_het[k*2]);
                i_to_j_d += __builtin_popcount(j_hom[k*2]);
            }

            if (j!=i+1)
                fprintf(f,"\t");
            fprintf(f,"%u", i_to_j_d);
        }
        fprintf(f,"\n");
    }

    fclose(f);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}
//}}}

//{{{ uint32_t wahbm_shared_by_name(char *in, char *out)
uint32_t wahbm_shared_by_name(char *in, char *out)
{

    FILE *f = fopen(out, "w");
    if (!f)
        err(EX_CANTCREAT, "Cannot open file \"%s\"", out);

    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_1_s, i_2_s, i_a_s, i_a_c_s;
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2 )
        err(EX_OSERR, "malloc error");

    uint32_t *i_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_a )
        err(EX_OSERR, "malloc error");
    uint32_t *i_a_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_a_c )
        err(EX_OSERR, "malloc error");

    uint32_t j_1_s, j_2_s, j_a_s;
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_a )
        err(EX_OSERR, "malloc error");

    uint32_t i,j,k, i_to_j_d;

    //for (i = 0; i < wf.num_records; ++i) {
    for (i = 0; i < wf->gqt_header->num_variants; ++i) {
        //load in the ith record
        //i_1_s = get_wah_bitmap_in_place(wf, i, 1, &i_1);
        i_1_s = get_wahbm_bitmap(wf, i, 1, &i_1);
        //i_2_s = get_wah_bitmap_in_place(wf, i, 2, &i_2);
        i_2_s = get_wahbm_bitmap(wf, i, 2, &i_2);

        memset(i_a, 0, sizeof(uint32_t)*max_wah_size);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_1, i_1_s);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_2, i_2_s);

        //for (j = i+1; j < wf.num_records; ++j) {
        for (j = i+1; j < wf->gqt_header->num_variants; ++j) {
            //load in the jth record
            //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
            j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
            //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
            j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);

            memset(j_a,   0, sizeof(uint32_t)*max_wah_size);

            j_a_s = wah_in_place_or(j_a, max_wah_size, j_1, j_1_s);
            j_a_s = wah_in_place_or(j_a, max_wah_size, j_2, j_2_s);

            j_a_s = wah_in_place_and(j_a, max_wah_size, i_a, i_a_s);

            i_to_j_d = 0;

            uint64_t *l_j_a = (uint64_t *)j_a;

            for (k = 0; k < j_a_s/2; k++) 
                i_to_j_d += __builtin_popcountll(l_j_a[k]);

            if (k*2 < j_a_s) 
                i_to_j_d += __builtin_popcount(j_a[k*2]);

            if (j!=i+1)
                fprintf(f,"\t");
            fprintf(f,"%u", i_to_j_d);
        }
        fprintf(f,"\n");
    }

    fclose(f);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}
//}}}

//{{{ uint32_t wahbm_shared_by_name(char *in, char *out)
uint32_t wahbm_shared_by_name_subpop(struct wahbm_file *wf,
                                     uint32_t *record_ids,
                                     uint32_t num_records)
{
    //struct wah_file wf = init_wahbm_file(in);

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_1_s, i_2_s, i_a_s, i_a_c_s;
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2 )
        err(EX_OSERR, "malloc error");

    uint32_t *i_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_a )
        err(EX_OSERR, "malloc error");
    uint32_t *i_a_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_a_c )
        err(EX_OSERR, "malloc error");

    uint32_t j_1_s, j_2_s, j_a_s;
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_a = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_a )
        err(EX_OSERR, "malloc error");

    uint32_t i,j,k, i_to_j_d;

    for (i = 0; i < num_records; ++i) {
        //load in the ith record
        //i_1_s = get_wah_bitmap_in_place(*wf, record_ids[i], 1, &i_1);
        i_1_s = get_wahbm_bitmap(wf, record_ids[i], 1, &i_1);
        //i_2_s = get_wah_bitmap_in_place(*wf, record_ids[i], 2, &i_2);
        i_2_s = get_wahbm_bitmap(wf, record_ids[i], 2, &i_2);

        memset(i_a, 0, sizeof(uint32_t)*max_wah_size);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_1, i_1_s);
        i_a_s = wah_in_place_or(i_a, max_wah_size, i_2, i_2_s);

        for (j = i+1; j < num_records; ++j) {
            //load in the jth record
            //j_1_s = get_wah_bitmap_in_place(*wf, record_ids[j], 1, &j_1);
            j_1_s = get_wahbm_bitmap(wf, record_ids[j], 1, &j_1);
            //j_2_s = get_wah_bitmap_in_place(*wf, record_ids[j], 2, &j_2);
            j_2_s = get_wahbm_bitmap(wf, record_ids[j], 2, &j_2);

            memset(j_a,   0, sizeof(uint32_t)*max_wah_size);

            j_a_s = wah_in_place_or(j_a, max_wah_size, j_1, j_1_s);
            j_a_s = wah_in_place_or(j_a, max_wah_size, j_2, j_2_s);

            j_a_s = wah_in_place_and(j_a, max_wah_size, i_a, i_a_s);

            i_to_j_d = 0;

            uint64_t *l_j_a = (uint64_t *)j_a;

            for (k = 0; k < j_a_s/2; k++) 
                i_to_j_d += __builtin_popcountll(l_j_a[k]);

            if (k*2 < j_a_s) 
                i_to_j_d += __builtin_popcount(j_a[k*2]);

            if (j!=i+1)
                printf("\t");
            printf("%u", i_to_j_d);
        }
        printf("\n");
    }

    return 0;
}
//}}}

//{{{ uint32_t wahbm_top_n_matches_by_name(char *in, uint32_t n)
uint32_t wahbm_top_n_matches_by_name(char *in, uint32_t num_matches)
{
    struct wahbm_file *wf = open_wahbm_file(in);

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;

    uint32_t i_0_s, i_1_s, i_2_s;
    uint32_t *i_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *i_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!i_2 )
        err(EX_OSERR, "malloc error");

    uint32_t j_0_s, j_1_s, j_2_s;
    uint32_t *j_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *j_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!j_2 )
        err(EX_OSERR, "malloc error");

    uint32_t x_0_s, x_1_s, x_2_s;
    uint32_t *x_1 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1 )
        err(EX_OSERR, "malloc error");
    uint32_t *x_2 = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2 )
        err(EX_OSERR, "malloc error");


    uint32_t x_0_sc, x_1_sc, x_2_sc;
    uint32_t *x_1_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_1_c )
        err(EX_OSERR, "malloc error");
    uint32_t *x_2_c = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!x_2_c )
        err(EX_OSERR, "malloc error");

    uint32_t t_s;
    uint32_t *t = (uint32_t *) malloc(sizeof(uint32_t)*max_wah_size);
    if (!t)
        err(EX_OSERR, "malloc error");

    uint32_t i,j,k,l, i_to_j_d;
    uint32_t leading_zeros, v;

    for (i = 0; i < wf->gqt_header->num_variants; ++i) {
        //load in the ith record

        //i_1_s = get_wah_bitmap_in_place(wf, i, 1, &i_1);
        i_1_s = get_wahbm_bitmap(wf, i, 1, &i_1);
        //i_2_s = get_wah_bitmap_in_place(wf, i, 2, &i_2);
        i_2_s = get_wahbm_bitmap(wf, i, 2, &i_2);

        memset(x_1_c, 0, sizeof(uint32_t)*max_wah_size);
        memset(x_2_c, 0, sizeof(uint32_t)*max_wah_size);

        x_1_sc = wah_in_place_or(x_1_c, max_wah_size, i_1, i_1_s);
        x_2_sc = wah_in_place_or(x_2_c, max_wah_size, i_2, i_2_s);
        
        for (j = i+1; j < wf->gqt_header->num_variants; ++j) {
            memcpy(x_1, x_1_c, x_1_sc * sizeof(uint32_t));
            memcpy(x_2, x_2_c, x_2_sc * sizeof(uint32_t));
            x_1_s = x_1_sc;
            x_2_s = x_2_sc;

            memset(j_1, 0, sizeof(uint32_t)*max_wah_size);
            memset(j_2, 0, sizeof(uint32_t)*max_wah_size);
            //load in the jth record
            //j_1_s = get_wah_bitmap_in_place(wf, j, 1, &j_1);
            j_1_s = get_wahbm_bitmap(wf, j, 1, &j_1);
            //j_2_s = get_wah_bitmap_in_place(wf, j, 2, &j_2);
            j_2_s = get_wahbm_bitmap(wf, j, 2, &j_2);
            //find the xor of all 4
            
            x_1_s =  wah_in_place_and(x_1, x_1_s, j_1, j_1_s);
            x_2_s =  wah_in_place_and(x_2, x_2_s, j_2, j_2_s);
            x_2_s =  wah_in_place_or(x_2, x_2_s, x_1, x_1_s);

            if (j != i+1)
                printf("\t");

            double score = 0.0;
            uint32_t num_hits = 0;
            for (k = 0; k < x_2_s; k++) {
#if 0
                if ( x_2[k] != 0 ) { 
                     score += __builtin_popcount(x_2[k]);
                }
#endif 
#if 1
                if ( x_2[k] != 0 ) { // this one has some number of hits
                    v = x_2[k];
                    //printf(" v%u:p%u:", v,__builtin_popcount(v));
                    for (l = 0; l < __builtin_popcount(v); ++l) {
                        leading_zeros = __builtin_clz(v);
                        //printf(" %u ", leading_zeros);
                        score = score + 
                            pow((((double)(k*31 + leading_zeros))/
                            ((double)wf->gqt_header->num_variants)),2);
                        //printf("L%u:", leading_zeros);
                        v &= ~(1 << (32 - leading_zeros - 1));
                        num_hits += 1;
                        if (num_hits == num_matches)
                            break;
                    }
                    if (num_hits == num_matches)
                        break;
                }
#endif
            }
#if 1
            if (num_hits == 0)
                printf("%f", 0.0);
            else
                printf("%f",score);
#endif
#if 0
                printf("%f",score);
#endif
        }
        printf("\n");
    }

    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);
    return 0;
}
//}}}

//{{{ uint32_t print_wahbm(struct wahbm_file *wf,
uint32_t print_wahbm(struct wahbm_file *wf,
                     uint32_t *record_ids,
                     uint32_t num_r,
                     uint32_t format)
{
    uint32_t i,j,k,l, bm_size, to_print = num_r;
    uint32_t *bm = NULL;

    uint32_t num_ints_per_record = 1 + ((wf->gqt_header->num_variants-1)/ 16);

    uint32_t *output = (uint32_t *)
            malloc((wf->gqt_header->num_variants)*sizeof(uint32_t));
    if (!output )
        err(EX_OSERR, "malloc error");

    if (num_r == 0)
        to_print = wf->gqt_header->num_samples;

    for (i = 0; i < to_print; ++i) {
        memset(output, 0, (wf->gqt_header->num_variants)*sizeof(uint32_t));
        for (j = 0; j < 4; ++j) {
            // get the compressed bitmap
            if (num_r > 0)
                bm_size = get_wahbm_bitmap(wf,
                                           record_ids[i],
                                           j,
                                           &bm);
            else
                bm_size = get_wahbm_bitmap(wf,
                                           i,
                                           j,
                                           &bm);

            // decompress 
            uint32_t *ints = NULL;
            uint32_t ints_size = wah_to_ints(bm,bm_size,&ints);
            uint32_t output_i = 0;
            free(bm);
            bm = NULL;

            for (k = 0; k < ints_size; ++k) {
                for (l = 0; l < 32; ++l) {
                    uint32_t bit = (ints[k] >> (31 - l)) & 1;
                    if (bit == 1) 
                        output[output_i] = j;

                    output_i += 1;

                    if (output_i >= wf->gqt_header->num_variants)
                        break;
                }
                if (output_i >= wf->gqt_header->num_variants)
                    break;
            }

            free(ints);
        }

        for (j = 0; j < wf->gqt_header->num_variants; ++j) {
            if (j != 0)
                printf(" ");
            printf("%u", output[j]);
        }
        printf("\n");

 

    }

#if 0
            // loop through each bit, and set the corresponding possition to j
            // if the bit is one
            int int_i = 0, bit_i = 0;
            for (k = 0; k < ints_size; ++k) {
                fprintf(stderr, "k:%u\n", k);
                for (l = 0; l < 32; ++l) {
                    uint32_t bit = (ints[k] >> (31 - l)) & 1;

                    if (bit == 1) {
                        tmp_record[int_i] += j << (30 - (bit_i * 2));
                    }

                    bit_i += 1;
                    if (bit_i == 16) {
                        int_i += 1;
                        bit_i = 0;
                    }

                    if (int_i >= num_ints_per_record)
                        break;
                }
                if (int_i >= num_ints_per_record)
                    break;
            }

            free(bm);
            free(ints);
            bm = NULL;
            ints = NULL;

            for (k = 0; k < num_ints_per_record; ++k) 
                output_record[k] += tmp_record[k];
        }

        uint32_t printed_bits = 0;
        for (j = 0; j < num_ints_per_record; ++j) {
            if (j !=0)
                printf(" ");
            for (k = 0; k < 16; ++k) {
                uint32_t bit = (output_record[j] >> (30 - 2*k)) & 3;
                if (k !=0)
                    printf(" ");
                printf("%u", bit);
                printed_bits += 1;
                if (printed_bits == wf.num_fields)
                    break;
            }
        }
        printf("\n");

    free(tmp_record);
    free(output_record);
#endif
    free(output);

    return to_print;
}
//}}}

//{{{ uint32_t print_by_name_wahbm(char *wahbm_file_name,
uint32_t print_by_name_wahbm(char *wahbm_file_name,
                               uint32_t *record_ids,
                               uint32_t num_r,
                               uint32_t format)
{
    //struct wah_file wf = init_wahbm_file(wahbm_file_name);
    struct wahbm_file *wf = open_wahbm_file(wahbm_file_name);
    uint32_t r = print_wahbm(wf, record_ids, num_r, format);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);
    return r;
}
//}}} 

#if 0
//{{{ uint32_t get_wah_bitmap(struct wahbm_file *wf,
uint32_t get_wah_bitmap(struct wah_file wf,
                        uint32_t wah_record,
                        uint32_t bitmap,
                        uint32_t **wah_bitmap)
{
    // get the size of the WAH-encoded bitmap
    uint64_t wah_size = 0;
    uint64_t wah_offset = 0;
    //fprintf(stderr, "wah_record:%u\tbitmap:%u\n", wah_record, bitmap);
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
                     sizeof(uint32_t) * 
                        (wf.record_offsets[wah_record*4 + bitmap] - wah_size);

        /*
        fprintf(stderr, "from:%llu\tto:%llu\tsize:%llu\n", 
                        wf.record_offsets[wah_record*4 + bitmap - 1],
                        wf.record_offsets[wah_record*4 + bitmap],
                        wah_size);
        */


    }
    //fprintf(stderr, "offset:%llu\n", wah_offset); 

    //fprintf(stderr, "wah_size:%llu\twah_offset:%llu\n", wah_size, wah_offset);


    *wah_bitmap = (uint32_t *) malloc(sizeof(uint32_t)*wah_size);
    if (!*wah_bitmap)
        err(EX_OSERR, "malloc error");
    fseek(wf.file, wah_offset, SEEK_SET);
    size_t fr = fread(*wah_bitmap,sizeof(uint32_t),wah_size,wf.file);
    check_file_read(wf.file_name, wf.file, wah_size, fr);

    return (uint32_t)wah_size;
}
//}}}
#endif

//{{{ uint32_t range_records_w_exclude_wahbm(struct wahbm_file *wf,
uint32_t range_records_w_exclude_wahbm(struct wahbm_file *wf,
                                           uint32_t *record_ids,
                                           uint32_t num_r,
                                           uint32_t start_test_value,
                                           uint32_t end_test_value,
                                           uint32_t exclude_value,
                                           uint32_t **R) 

{
    uint32_t *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    uint32_t record_curr_bm_size = 0,
                 record_new_bm_size = 0,
                 record_tmp_bm_size = 0;

    uint32_t *query_curr_bm = NULL,
                 *query_tmp_bm = NULL;

    uint32_t query_curr_bm_size = 0,
                 query_tmp_bm_size = 0;


    uint32_t i,j,k,l;

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

            record_new_bm_size = get_wahbm_bitmap(wf,
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

//{{{ uint32_t count_range_records_wahbm(struct wahbm_file *wf,
uint32_t count_range_records_wahbm(struct wahbm_file *wf,
                                       uint32_t *record_ids,
                                       uint32_t num_r,
                                       uint32_t start_test_value,
                                       uint32_t end_test_value,
                                       uint32_t **R) 

{
    //*R = (uint32_t *) calloc(wf.num_fields,sizeof(uint32_t));
    *R = (uint32_t *) calloc(wf->gqt_header->num_variants,sizeof(uint32_t));

    uint32_t *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    uint32_t record_curr_bm_size = 0,
                 record_new_bm_size = 0,
                 record_tmp_bm_size = 0;

    uint32_t i,j, r_size;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);
            */
            record_new_bm_size = get_wahbm_bitmap(wf,
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
                           wf->gqt_header->num_variants,
                           record_curr_bm,
                           record_curr_bm_size);
    }
    return wf->gqt_header->num_variants;
}
//}}}

//{{{ uint32_t sum_range_records_wahbm(struct wahbm_file *wf,
uint32_t sum_range_records_wahbm(struct wahbm_file *wf,
                                 uint32_t *record_ids,
                                 uint32_t num_r,
                                 uint32_t **R) 

{
    *R = (uint32_t *) calloc(wf->gqt_header->num_variants,sizeof(uint32_t));

    uint32_t *record_curr_bm = NULL,
             *record_new_bm = NULL,
             *record_tmp_bm = NULL;

    uint32_t record_curr_bm_size = 0,
             record_new_bm_size = 0,
             record_tmp_bm_size = 0;

    uint32_t i,j, r_size;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = 0; j < 4; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);
            */
            record_new_bm_size = get_wahbm_bitmap(wf,
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
                           wf->gqt_header->num_variants,
                           record_curr_bm,
                           record_curr_bm_size);
    }
    return wf->gqt_header->num_variants;
}
//}}}

//{{{ uint32_t range_records_wahbm(struct wahbm_file *wf,
uint32_t range_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t start_test_value,
                              uint32_t end_test_value,
                              uint32_t **R) 

{
    uint32_t *record_curr_bm = NULL,
                 *record_new_bm = NULL,
                 *record_tmp_bm = NULL;

    uint32_t record_curr_bm_size = 0,
                 record_new_bm_size = 0,
                 record_tmp_bm_size = 0;

    uint32_t *query_curr_bm = NULL,
                 *query_tmp_bm = NULL;

    uint32_t query_curr_bm_size = 0,
                 query_tmp_bm_size = 0;


    uint32_t i,j,k,l;

    for (i = 0; i < num_r; ++i) {
        // or all of the bit maps for this record then and that will a 
        // running total

        record_curr_bm = NULL;
        record_new_bm = NULL;
        record_tmp_bm = NULL;

        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap(wf,
                                                record_ids[i],
                                                j,
                                                &record_new_bm);
            */
            record_new_bm_size = get_wahbm_bitmap(wf,
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

//{{{ uint32_t add_wahbm(uint32_t *R,
uint32_t add_wahbm(uint32_t *R,
                   uint32_t r_size,
                   uint32_t *wah,
                   uint32_t wah_size)
{

    uint32_t wah_c,
             wah_i,
             num_words,
             fill_bit,
             bits,
             bit,
             bit_i,
             word_i,
             field_i;

    field_i = 0;

    uint32_t v;

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
            field_i += num_words * 31;
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
//{{{void avx_add(uint32_t bits,
void avx_add(uint32_t bits,
             __m256i *s_1,
             __m256i *s_2,
             __m256i *s_3,
             __m256i *s_4,
             __m256i *m,
             __m256i *R_avx,
             uint32_t field_i)
{

    uint32_t avx_i = field_i/8;

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
//{{{ uint32_t avx_add_wahbm(uint32_t *R,
uint32_t avx_add_wahbm(uint32_t *R,
                       uint32_t r_size,
                       uint32_t *wah,
                       uint32_t wah_size)
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

    uint32_t wah_c,
                 wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    field_i = 0;

    uint32_t buf, buf_empty_bits = 32;

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
                            uint32_t padding = field_i % 32;
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

#ifdef __AVX2__
//{{{void avx_add_n(uint32_t bits,
void avx_add_n(uint32_t bits,
             __m256i *s_1,
             __m256i *s_2,
             __m256i *s_3,
             __m256i *s_4,
             __m256i *m,
             __m256i *N,
             __m256i *R_avx,
             uint32_t field_i)
{
    uint32_t avx_i = field_i/8;

    __m256i y1 = _mm256_set1_epi32(bits);

    __m256i y2 = _mm256_srlv_epi32 (y1, *s_1);
    __m256i y3 = _mm256_and_si256 (y2, *m);
    __m256i y4 = _mm256_mullo_epi16(y3, *N);
    R_avx[3+avx_i] = _mm256_add_epi32(R_avx[3+avx_i], y4);

    y2 = _mm256_srlv_epi32 (y1, *s_2);
    y3 = _mm256_and_si256 (y2, *m);
    y4 = _mm256_mullo_epi16(y3, *N);
    R_avx[2+avx_i] = _mm256_add_epi32(R_avx[2+avx_i], y4);

    y2 = _mm256_srlv_epi32 (y1, *s_3);
    y3 = _mm256_and_si256 (y2, *m);
    y4 = _mm256_mullo_epi16(y3, *N);
    R_avx[1+avx_i] = _mm256_add_epi32(R_avx[1+avx_i], y4);

    y2 = _mm256_srlv_epi32 (y1, *s_4);
    y3 = _mm256_and_si256 (y2, *m);
    y4 = _mm256_mullo_epi16(y3, *N);
    R_avx[0+avx_i] = _mm256_add_epi32(R_avx[0+avx_i], y4);
}
//}}}
#endif

#ifdef __AVX2__
//{{{ uint32_t avx_add_n_wahbm(uint32_t *R,
uint32_t avx_add_n_wahbm(uint32_t *R,
                             uint32_t n,
                             uint32_t r_size,
                             uint32_t *wah,
                             uint32_t wah_size)
{
    __attribute__((aligned(64))) int rshift_4[8] = 
            { 31, 30, 29, 28, 27, 26, 25, 24 };
    __attribute__((aligned(64))) int rshift_3[8] =
            { 23, 22, 21, 20, 19, 18, 17, 16 };
    __attribute__((aligned(64))) int rshift_2[8] =
            { 15, 14, 13, 12, 11, 10, 9, 8 };
    __attribute__((aligned(64))) int rshift_1[8] =
            { 7, 6, 5, 4, 3, 2, 1, 0 };
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

    __m256i n_avx = _mm256_set1_epi32(n);

    uint32_t wah_c,
                 wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    field_i = 0;

    uint32_t buf, buf_empty_bits = 32;

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
                avx_add_n(buf,
                          &s_1,
                          &s_2,
                          &s_3,
                          &s_4,
                          &m,
                          &n_avx,
                          R_avx,
                          field_i);
            
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
                    avx_add_n(buf,
                              &s_1,
                              &s_2,
                              &s_3,
                              &s_4,
                              &m,
                              &n_avx,
                              R_avx,
                              field_i);
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
                            uint32_t padding = field_i % 32;
                            buf = bits >> (padding - 1);
                            avx_add_n(buf,
                                      &s_1,
                                      &s_2,
                                      &s_3,
                                      &s_4,
                                      &m,
                                      &n_avx,
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

                        avx_add_n(buf,
                                  &s_1,
                                  &s_2,
                                  &s_3,
                                  &s_4,
                                  &m,
                                  &n_avx,
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
        R[field_i] += n*((buf >> (31 - bit_i)) & 1);
        field_i += 1;

        if (field_i >= r_size)
            return r_size;
    }

    return r_size;
}
//}}}
#endif

//{{{ uint32_t add_n_wahbm(uint32_t *R,
uint32_t add_n_wahbm(uint32_t *R,
                       uint32_t n,
                       uint32_t r_size,
                       uint32_t *wah,
                       uint32_t wah_size)
{

    uint32_t wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    uint32_t t;
    field_i = 0;

    for (wah_i = 0; wah_i < wah_size; ++wah_i) {
        /* From the current  value:
         * 1) determine how many works (1 if it is a litteral and more than one
         * if it is a fill.
         * 2) get the bits to add, if it is a fill then take the fill bit and
         * create a word of only that bit, other wiswe grab the literal
         */
        if (wah[wah_i] >> 31 == 1) {
            num_words = (wah[wah_i] & 0x3fffffff);
            fill_bit = (wah[wah_i]>=0xC0000000?1:0);
            bits = (fill_bit?0x7FFFFFFF:0);
        } else {
            num_words = 1;
            bits = wah[wah_i];
        }

        // If there is nothing to add for more than one word, skip
        if ( (num_words > 1) && (fill_bit == 0) ) {
            field_i += num_words * 31;
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

    uint32_t bit_i,
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

    uint32_t bit_i,
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

//{{{ uint32_t p_pool_add_n_wahbm(uint32_t *R,
uint32_t p_pool_add_n_wahbm(uint32_t *R,
                                uint32_t n,
                                uint32_t r_size,
                                uint32_t *wah,
                                uint32_t wah_size,
                                struct pool *t_pool)
{
    uint32_t wah_i,
                 num_words,
                 fill_bit,
                 bits,
                 bit,
                 bit_i,
                 word_i,
                 field_i;
    uint32_t t;
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
            field_i += num_words * 31;
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
                    if (!arg)
                        err(EX_OSERR, "malloc error");
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
//{{{ uint32_t eq_records_wahbm(struct wahbm_file *wf,
uint32_t eq_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value, test_value+1, R);
}
//}}}

//{{{ uint32_t ne_records_wahbm(struct wahbm_file *wf,
uint32_t ne_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for lower bound and upper bound.
    // exclude the test_value
    return range_records_w_exclude_wahbm(wf, record_ids, num_r, 0, 4, test_value, R);
}
//}}}

//{{{ uint32_t gt_records_wahbm(struct wahbm_file *wf,
uint32_t gt_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value+1, 4, R);
}
//}}}

//{{{ uint32_t gte_records_wahbm(struct wahbm_file *wf,
uint32_t gte_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, test_value, 4, R);
}
//}}}

//{{{ uint32_t lt_records_wahbm(struct wahbm_file *wf,
uint32_t lt_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, 0, test_value, R);
}
//}}}

//{{{ uint32_t lte_records_wahbm(struct wahbm_file *wf,
uint32_t lte_records_wahbm(struct wahbm_file *wf,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t test_value,
                              uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return range_records_wahbm(wf, record_ids, num_r, 0, test_value+1, R);
}
//}}}
//}}}

//{{{ uint32_t gt_count_records_wahbm(struct wahbm_file *wf,
uint32_t gt_count_records_wahbm(struct wahbm_file *wf,
                                    uint32_t *record_ids,
                                    uint32_t num_r,
                                    uint32_t test_value,
                                    uint32_t **R)
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
