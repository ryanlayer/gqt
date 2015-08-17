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

#include "wahbm_in_place.h"
#include "wah.h"
#include "wahbm.h"

#include "timer.h"
#include "pthread_pool.h"


// wahbm in place
//{{{ uint32_t wah_in_place_or(uint32_t *r_wah,
uint32_t wah_in_place_or(uint32_t *r_wah,
                         uint32_t r_wah_size,
                         uint32_t *wah,
                         uint32_t wah_size)
{

    uint32_t r_wah_i = 0, wah_c = 0;
    uint32_t wah_i, fill_size, wah_v, end;
    for (wah_i = 0; wah_i < wah_size; ++wah_i)
    {
        wah_v = wah[wah_i];
        // is the current word a fill
        if (wah_v >= 0x80000000) {
            fill_size = wah_v & 0x3fffffff;
            wah_c += fill_size;
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
            wah_c += 1;
        }
    }

    return r_wah_size;
}
//}}}

//{{{ uint32_t wah_in_place_and(uint32_t *r_wah,
uint32_t wah_in_place_and(uint32_t *r_wah,
                               uint32_t r_wah_size,
                               uint32_t *wah,
                               uint32_t wah_size)
{

    uint32_t r_wah_i = 0;
    uint32_t wah_i, fill_size, wah_v, end;
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

//{{{ uint32_t wah_in_place_xor(uint32_t *r_wah,
uint32_t wah_in_place_xor(uint32_t *r_wah,
                          uint32_t r_wah_size,
                          uint32_t *wah,
                          uint32_t wah_size)
{

    uint32_t r_wah_i = 0;
    uint32_t wah_i, fill_size, fill_v, wah_v, end;
    for (wah_i = 0; wah_i < wah_size; ++wah_i)
    {
        wah_v = wah[wah_i];
        // is the current word a fill
        if (wah_v >= 0x80000000) {
            fill_size = wah_v & 0x3fffffff;
            if (wah_v >> 30 == 3) // fill of 1s
                fill_v = 0x7fffffff;
            else  // fill of 0s
                fill_v = 0;

            end =  r_wah_i + fill_size;
            for ( ; r_wah_i < end; ++r_wah_i)
                r_wah[r_wah_i] = r_wah[r_wah_i] ^ fill_v;
        } else {
            r_wah[r_wah_i] = r_wah[r_wah_i] ^ wah[wah_i];
            r_wah_i += 1;
        }
    }

    return r_wah_i;
}
//}}}

//{{{ uint32_t range_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t range_records_in_place_wahbm(struct wahbm_file *wf,
                                      uint32_t *record_ids,
                                      uint32_t num_r,
                                      uint32_t start_test_value,
                                      uint32_t end_test_value,
                                      uint32_t **R) 
{

    //uint32_t max_wah_size = (wf.num_fields + 31 - 1)/ 31;
    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
    uint32_t *record_new_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!record_new_bm )
        err(EX_OSERR, "malloc error");

    uint32_t *or_result_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!or_result_bm )
        err(EX_OSERR, "malloc error");
    uint32_t *and_result_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!and_result_bm )
        err(EX_OSERR, "malloc error");
    uint32_t and_result_bm_size, record_new_bm_size, or_result_bm_size;
    uint32_t i,j;
    for (i = 0; i < max_wah_size; ++i)
        and_result_bm[i] = 0x7fffffff;

    for (i = 0; i < num_r; ++i) {
        // or the appropriate bitmaps
        memset(or_result_bm, 0, sizeof(uint32_t)*max_wah_size);

        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);
            */

            record_new_bm_size = get_wahbm_bitmap(wf,
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

//{{{ uint32_t count_range_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t count_range_records_in_place_wahbm(struct wahbm_file *wf,
                                            uint32_t *record_ids,
                                            uint32_t num_r,
                                            uint32_t start_test_value,
                                            uint32_t end_test_value,
                                            uint32_t **R) 
{
    //*R = (uint32_t *) calloc(wf.num_fields,sizeof(uint32_t));
    *R = (uint32_t *) calloc(wf->gqt_header->num_variants,sizeof(uint32_t));
    if (!*R )
        err(EX_OSERR, "malloc error");

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
    uint32_t *record_new_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!record_new_bm )
        err(EX_OSERR, "malloc error");

    uint32_t *or_result_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!or_result_bm )
        err(EX_OSERR, "malloc error");

    uint32_t and_result_bm_size, record_new_bm_size, or_result_bm_size;
    uint32_t i,j,r_size;

    for (i = 0; i < num_r; ++i) {
        // or the appropriate bitmaps
        memset(or_result_bm, 0, sizeof(uint32_t)*max_wah_size);

        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);
            */
            record_new_bm_size = get_wahbm_bitmap(wf,
                                                  record_ids[i],
                                                  j,
                                                  &record_new_bm);


            or_result_bm_size = wah_in_place_or(or_result_bm,
                                                max_wah_size,
                                                record_new_bm,
                                                record_new_bm_size); 
        }

        r_size = add_wahbm(*R,
                           wf->gqt_header->num_variants,
                           or_result_bm,
                           or_result_bm_size);
    }

    free(record_new_bm);
    free(or_result_bm);

    return wf->gqt_header->num_variants;
}
//}}}

#ifdef __AVX2__
//{{{ uint32_t avx_count_range_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t avx_count_range_records_in_place_wahbm(struct wahbm_file *wf,
                                                uint32_t *record_ids,
                                                uint32_t num_r,
                                                uint32_t start_test_value,
                                                uint32_t end_test_value,
                                                uint32_t **R) 

{
    //*R = (uint32_t *) calloc(wf.num_fields,sizeof(uint32_t));
    //int r = posix_memalign((void **)R, 32, wf.num_fields*sizeof(uint32_t));
    int r = posix_memalign((void **)R,
                            32,
                            wf->gqt_header->num_variants*sizeof(uint32_t));
    if (r !=0 )
        err(EX_OSERR, "posix_memalign error");
    //memset(*R, 0, wf.num_fields*sizeof(uint32_t));
    memset(*R, 0, (wf->gqt_header->num_variants)*sizeof(uint32_t));

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
    uint32_t *record_new_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!record_new_bm )
        err(EX_OSERR, "malloc error");

    uint32_t *or_result_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!or_result_bm )
        err(EX_OSERR, "malloc error");

    uint32_t and_result_bm_size, record_new_bm_size, or_result_bm_size;
    uint32_t i,j,r_size;

    for (i = 0; i < num_r; ++i) {
        // or the appropriate bitmaps
        memset(or_result_bm, 0, sizeof(uint32_t)*max_wah_size);

        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);
            */
            record_new_bm_size = get_wahbm_bitmap(wf,
                                                  record_ids[i],
                                                  j,
                                                  &record_new_bm);

            or_result_bm_size = wah_in_place_or(or_result_bm,
                                                max_wah_size,
                                                record_new_bm,
                                                record_new_bm_size); 
        }

       r_size = avx_add_wahbm(*R,
                               wf->gqt_header->num_variants,
                               or_result_bm,
                               or_result_bm_size);
    }

    free(record_new_bm);
    free(or_result_bm);

    return wf->gqt_header->num_variants;
}
//}}}
#endif

#ifdef __AVX2__
//{{{ uint32_t avx_sum_range_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t avx_sum_range_records_in_place_wahbm(struct wahbm_file *wf,
                                              uint32_t *record_ids,
                                              uint32_t num_r,
                                              uint32_t start_test_value,
                                              uint32_t end_test_value,
                                              uint32_t **R)

{
    //int r = posix_memalign((void **)R, 32, wf.num_fields*sizeof(uint32_t));
    int r = posix_memalign((void **)R,
                           32,
                           (wf->gqt_header->num_variants)*sizeof(uint32_t));
    if (r !=0 )
        err(EX_OSERR, "posix_memalign error");
    //memset(*R, 0, wf.num_fields*sizeof(uint32_t));
    memset(*R, 0, (wf->gqt_header->num_variants)*sizeof(uint32_t));

    //uint32_t max_wah_size = (wf.num_fields + 31 - 1)/ 31;
    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
    uint32_t *record_new_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!record_new_bm )
        err(EX_OSERR, "malloc error");

    uint32_t record_new_bm_size;
    uint32_t i,j,r_size;
#ifdef time_sum_range_records_in_place_wahbm
    unsigned long t1 = 0, t2 = 0, t3 = 0;
#endif

    for (i = 0; i < num_r; ++i) {
        for (j = start_test_value; j < end_test_value; ++j) {

#ifdef time_sum_range_records_in_place_wahbm
            start();
#endif
            /*
            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);
            */

            record_new_bm_size = get_wahbm_bitmap(wf,
                                                  record_ids[i],
                                                  j,
                                                  &record_new_bm);

#ifdef time_sum_range_records_in_place_wahbm
            stop();
            t1+=report();
#endif

            r_size = avx_add_n_wahbm(*R,
                                     j,
                                     //wf.num_fields,
                                     wf->gqt_header->num_variants,
                                     record_new_bm,
                                     record_new_bm_size);

#ifdef time_sum_range_records_in_place_wahbm
            stop();
            t2+=report();
#endif
        }
    }

#ifdef time_sum_range_records_in_place_wahbm
    unsigned long tall = t1 + t2;
    fprintf(stderr,"%lu %f\t"
                   "%lu %f\t"
                   "%lu\n", 
            t1,
            ((double)t1)/((double)tall),
            t2,
            ((double)t2)/((double)tall),
            tall);

#endif
    free(record_new_bm);
    //return wf.num_fields;
    return wf->gqt_header->num_variants;
}
//}}}
#endif

//{{{ uint32_t sum_range_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t sum_range_records_in_place_wahbm(struct wahbm_file *wf,
                                          uint32_t *record_ids,
                                          uint32_t num_r,
                                          uint32_t start_test_value,
                                          uint32_t end_test_value,
                                          uint32_t **R)

{
    //*R = (uint32_t *) calloc(wf.num_fields,sizeof(uint32_t));
    *R = (uint32_t *) calloc(wf->gqt_header->num_variants,sizeof(uint32_t));
    if (!*R )
        err(EX_OSERR, "malloc error");

    uint32_t max_wah_size = (wf->gqt_header->num_variants + 31 - 1)/ 31;
    uint32_t *record_new_bm = (uint32_t *)
                        malloc(sizeof(uint32_t)*max_wah_size);
    if (!record_new_bm )
        err(EX_OSERR, "malloc error");

    uint32_t and_result_bm_size, record_new_bm_size, or_result_bm_size;
    uint32_t i,j,r_size;

    struct pool *t_pool = pool_start(t_add_n_wahbm_2, 2);

    for (i = 0; i < num_r; ++i) {
        for (j = start_test_value; j < end_test_value; ++j) {

            /*
            record_new_bm_size = get_wah_bitmap_in_place(wf,
                                                         record_ids[i],
                                                         j,
                                                         &record_new_bm);
            */

            record_new_bm_size = get_wahbm_bitmap(wf,
                                                  record_ids[i],
                                                  j,
                                                  &record_new_bm);

            r_size = add_n_wahbm(*R,
                                 j,
                                 wf->gqt_header->num_variants,
                                 record_new_bm,
                                 record_new_bm_size);

#if 0
            r_size = p_pool_add_n_wahbm(*R,
                                        j,
                                        wf.num_fields,
                                        record_new_bm,
                                        record_new_bm_size,
                                        t_pool);
            pool_wait(t_pool);
#endif
        }

    }
    free(record_new_bm);
    return wf->gqt_header->num_variants;
}
//}}}

//{{{ uint32_t gt_count_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t gt_count_records_in_place_wahbm(struct wahbm_file *wf,
                                         uint32_t *record_ids,
                                         uint32_t num_r,
                                         uint32_t test_value,
                                         uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return count_range_records_in_place_wahbm(wf,
                                              record_ids,
                                              num_r,
                                              test_value+1,
                                              4,
                                              R);
}
//}}}

#ifdef __AVX2__
//{{{ uint32_t avx_gt_count_records_in_place_wahbm(struct wahbm_file *wf,
uint32_t avx_gt_count_records_in_place_wahbm(struct wahbm_file *wf,
                                             uint32_t *record_ids,
                                             uint32_t num_r,
                                             uint32_t test_value,
                                             uint32_t **R) 

{
    // TODO: need constants for upper bound.
    return avx_count_range_records_in_place_wahbm(wf,
                                                  record_ids,
                                                  num_r,
                                                  test_value+1,
                                                  4,
                                                  R);
}
//}}}
#endif

#if 0
//{{{ uint32_t get_wah_bitmaps_in_place(struct wah_file wf,
uint32_t get_wah_bitmaps_in_place(struct wah_file wf,
                                      uint32_t wah_record,
                                      uint32_t **wah_bitmap,
                                      uint32_t *wah_sizes)
{
    // get the size of the WAH-encoded bitmap
    uint64_t wah_size = 0, wah_offset = 0;
    if (wah_record == 0) {
        wah_size = wf.record_offsets[wah_record + 3];
        wah_offset = wf.header_offset;
    } else {

        wah_size = wf.record_offsets[wah_record*4 + 3] - 
                   wf.record_offsets[wah_record*4 - 1];

        wah_offset = wf.header_offset +
                     sizeof(uint32_t) * 
                        (wf.record_offsets[wah_record*4 - 1]);
    }

    uint32_t i;
    for (i = 0; i < 4; ++i) {
        if ((wah_record == 0) && (i == 0)) 
            wah_sizes[i] = wf.record_offsets[wah_record];
        else
            wah_sizes[i] = wf.record_offsets[wah_record*4 + i] - 
                   wf.record_offsets[wah_record*4 + i - 1];
    }


    fseek(wf.file, wah_offset, SEEK_SET);
    size_t fr = fread(*wah_bitmap,sizeof(uint32_t),wah_size,wf.file);
    check_file_read(wf.file_name, wf.file, wah_size, fr);

    return (uint32_t)wah_size;
}
//}}}

//{{{ uint32_t get_wah_bitmap_in_place(struct wah_file wf,
uint32_t get_wah_bitmap_in_place(struct wah_file wf,
                                     uint32_t wah_record,
                                     uint32_t bitmap,
                                     uint32_t **wah_bitmap)
{
    // get the size of the WAH-encoded bitmap
    uint64_t wah_size = 0;
    uint64_t wah_offset = 0;
    if ((wah_record == 0) && (bitmap == 0)) {
        wah_size = wf.record_offsets[wah_record + bitmap];
        wah_offset = wf.header_offset;
    } else {
        wah_size = wf.record_offsets[wah_record*4 + bitmap] - 
                   wf.record_offsets[wah_record*4 + bitmap - 1];

        wah_offset = wf.header_offset +
                     sizeof(uint32_t) * 
                        (wf.record_offsets[wah_record*4 + bitmap] - wah_size);
    }

    fseek(wf.file, wah_offset, SEEK_SET);
    size_t fr = fread(*wah_bitmap,sizeof(uint32_t),wah_size,wf.file);
    check_file_read(wf.file_name, wf.file, wah_size, fr);

    return (uint32_t)wah_size;
}
//}}}
#endif
