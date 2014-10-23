#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include "genotq.h"
#include "timer.h"

// From http://www.hackersdelight.org/hdcodetxt/nlz.c.txt
//{{{int nlz1(unsigned x)
int nlz1(unsigned x)
{
    int n;

    if (x == 0) return(32);
    n = 0;
    if (x <= 0x0000FFFF) {n = n +16; x = x <<16;}
    if (x <= 0x00FFFFFF) {n = n + 8; x = x << 8;}
    if (x <= 0x0FFFFFFF) {n = n + 4; x = x << 4;}
    if (x <= 0x3FFFFFFF) {n = n + 2; x = x << 2;}
    if (x <= 0x7FFFFFFF) {n = n + 1;}
    return n;
}
//}}}

//{{{struct bcf_file init_bcf_file(char *file_name)
struct bcf_file init_bcf_file(char *file_name)
{
    struct bcf_file bcf_f;

    bcf_f.fp = hts_open(file_name,"rb");
    bcf_f.hdr = bcf_hdr_read(bcf_f.fp);
    bcf_f.line = bcf_init1();
    bcf_f.num_records = bcf_hdr_nsamples(bcf_f.hdr);
    bcf_f.gt = NULL;

    return bcf_f;
}
//}}}

//{{{int read_unpack_next_bcf_line(struct bcf_file *bcf_f,
int read_unpack_next_bcf_line(struct bcf_file *bcf_f,
                              int *num_samples,
                              int *num_gts_per_sample)
{
    int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);

    if (r < 0)
        return r; 

    int ntmp = bcf_f->line->n_sample;
    bcf_unpack(bcf_f->line, BCF_UN_ALL);
    *num_gts_per_sample = bcf_get_genotypes(bcf_f->hdr,
                                            bcf_f->line,
                                            &(bcf_f->gt),
                                            &ntmp);
    *num_samples = bcf_hdr_nsamples(bcf_f->hdr);
    *num_gts_per_sample /= *num_samples;

    return r;
}
//}}}

//{{{uint32_t pack_sum_count_prefix_bcf_line(struct bcf_file bcf_f,
uint32_t pack_sum_count_prefix_bcf_line(struct bcf_file bcf_f,
                                        uint32_t num_samples,
                                        uint32_t num_gts_per_sample,
                                        uint32_t **packed_ints,
                                        uint32_t *sum,
                                        uint32_t *prefix_len)
{
    uint32_t num_ints = 1 + ((num_samples - 1) / 16);
    *packed_ints = calloc(num_ints, sizeof(uint32_t));
    //fprintf(stderr, "pack_sum_count_prefix_bcf_line\tnum_ints:%u\n",num_ints); 
    int32_t *gt_i = bcf_f.gt;

    uint32_t two_bit_i = 0, int_i = 0;

    *sum = 0;

    uint32_t i, j, a = 0;
    for (i = 0; i < num_samples; ++i) {
        uint32_t gt = 0;
        for (j=0; j< num_gts_per_sample; ++j) {
            gt += bcf_gt_allele(gt_i[j]);
        }

        //fprintf(stderr, "pack_sum_count_prefix_bcf_line\tint_i:%u\n",int_i); 
        (*packed_ints)[int_i] += gt << (30 - 2*two_bit_i);
        two_bit_i += 1;
        if (two_bit_i == 16) {
            two_bit_i = 0;
            int_i += 1;
        }

        *sum += gt;
        gt_i += num_gts_per_sample;
    }

    *prefix_len = 0;
    for (i = 0; i < num_ints; ++i) {
        if ( (*packed_ints)[i] == 0 )
            *prefix_len += 32;
        else {
            *prefix_len += nlz1((*packed_ints)[i]);
            break;
        }
    }

    return num_ints;
}
//}}}

//{{{uint32_t md_bcf_line(struct bcf_file bcf_f,
uint32_t md_bcf_line(struct bcf_file bcf_f,
                     char **md)
{
    size_t len = strlen(bcf_hdr_id2name(bcf_f.hdr, bcf_f.line->rid)) +
                 10 + // max length of pos
                 strlen(bcf_f.line->d.id) +
                 strlen(bcf_f.line->d.allele[0]) +
                 strlen(bcf_f.line->d.allele[1]) +
                 4; //tabs
    *md = (char *) malloc(len * sizeof(char));

    sprintf(*md,
            "%s\t%d\t%s\t%s\t%s",
            bcf_hdr_id2name(bcf_f.hdr, bcf_f.line->rid),
            bcf_f.line->pos,
            bcf_f.line->d.id,
            bcf_f.line->d.allele[0],
            bcf_f.line->d.allele[1]);

    return len;
}
//}}}}

//{{{ void push_bcf_gt_md(pri_queue *q,
void push_bcf_gt_md(pri_queue *q,
                    struct bcf_file *bcf_f,
                    struct hdf5_file *hdf5_f)
{
    int num_samples, num_gts_per_sample;
    uint32_t *packed_ints;
    uint32_t sum, prefix_len, packed_ints_len;
    int r;
    char *md;
    uint32_t md_len;

    int i,j;

    priority p;

    for (i = 0; i < hdf5_f->num_vars; ++i) {

        r = read_unpack_next_bcf_line(bcf_f,
                                      &num_samples,
                                      &num_gts_per_sample);

        packed_ints_len = pack_sum_count_prefix_bcf_line(*bcf_f,
                                                         num_samples,
                                                         num_gts_per_sample,
                                                         &packed_ints,
                                                         &sum,
                                                         &prefix_len);
        md_len = md_bcf_line(*bcf_f, &md);

        write_hdf5_gt(*hdf5_f, i, packed_ints, md);

        p.sum = sum;
        p.len = prefix_len;
        int *j = (int *) malloc (sizeof(int));
        j[0] = i;
        priq_push(*q, j, p);

        free(packed_ints);
        free(md);
    }
}
//}}}

//{{{int convert_file_by_name_bcf_to_wahbm_bim(char *in,
int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out)
{

    struct hdf5_file hdf5_f = init_hdf5_file(".tmp.h5",
                                             num_records, //num_vars,
                                             num_fields); //num_inds);

    struct bcf_file bcf_f = init_bcf_file(in);

    pri_queue q = priq_new(0);

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    unsigned long t1 = 0, t2 = 0, t3 = 0;
#endif


#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    start();
#endif

    push_bcf_gt_md(&q, &bcf_f, &hdf5_f);

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    stop();
    t1 = report();
#endif

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    start();
#endif

    sort_rotate_gt_md(&q, &hdf5_f, bim_out);

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    stop();
    t2 = report();
#endif

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    start();
#endif

    int r = convert_hdf5_ind_ubin_to_ind_wah(hdf5_f, wah_out);

#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    stop();
    t3 = report();
#endif

    
#ifdef time_convert_file_by_name_bcf_to_wahbm_bim
    double t = t1 + t2 + t3 + 0.0;
    fprintf(stderr, "push_bcf_gt_md: %lu %f\t"
                    "sort_rotate_gt_md: %lu %f\t"
                    "convert_hdf5_ind_ubin_to_ind_wah: %lu %f\n",
                    t1, t1/t,
                    t2, t2/t,
                    t3, t3/t
           );
#endif

    priq_free(q);
    close_hdf5_file(hdf5_f);
    close_bcf_file(&bcf_f);
    remove(".tmp.h5");

    //return r;
    return 0;
}
//}}}

//{{{void close_bcf_file(struct bcf_file *bcf_f)
void close_bcf_file(struct bcf_file *bcf_f)
{
    bcf_hdr_destroy(bcf_f->hdr);
    bcf_destroy1(bcf_f->line);
    hts_close(bcf_f->fp);
}
//}}}
