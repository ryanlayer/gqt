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


//{{{int convert_file_by_name_bcf_to_wahbm_bim(char *in,
int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out)
{
    uint32_t num_inds = num_fields;
    uint32_t num_vars = num_records;
    char *gt_of_name = ".gt.tmp.packed";
    char *s_gt_of_name = ".s.gt.tmp.packed";
    char *r_s_gt_of_name = ".r.s.gt.tmp.packed";
    char *md_of_name = ".md.tmp.packed";

    struct bcf_file bcf_f = init_bcf_file(in);
    pri_queue q = priq_new(0);

    uint32_t *md_index = (uint32_t *) malloc(num_vars * sizeof(uint32_t));

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    sort_gt_md(&q,
               md_index,
               num_inds,
               num_vars,
               gt_of_name,
               s_gt_of_name,
               md_of_name,
               bim_out);

    rotate_encode_wahbm(num_inds,
                        num_vars,
                        s_gt_of_name,
                        r_s_gt_of_name);

    close_bcf_file(&bcf_f);

    return 0;
}
//}}}

//{{{ void push_bcf_gt_md(pri_queue *q,
void push_bcf_gt_md(pri_queue *q,
                    struct bcf_file *bcf_f,
                    uint32_t *md_index,
                    uint32_t num_inds,
                    uint32_t num_vars,
                    char *gt_of_name,
                    char *md_of_name)
{

    FILE *gt_of = fopen(gt_of_name,"wb");
    FILE *md_of = fopen(md_of_name,"w");

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t *packed_ints = (uint32_t *) calloc(num_ind_ints,
                                                sizeof(uint32_t));

    uint32_t md_i = 0;

    uint32_t i, j, k, sum, int_i, two_bit_i = 0;
    int ntmp = 0;

    int32_t *gt_p = NULL;

    priority p;

    for (i = 0; i < num_vars; ++i) {
        sum = 0;
        int_i = 0;
        two_bit_i = 0;

        // Get the next bcf record
        int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);

        // Unpack all of the fields 
        bcf_unpack(bcf_f->line, BCF_UN_ALL);

        // Get metadata
        size_t len = strlen(bcf_hdr_id2name(bcf_f->hdr, bcf_f->line->rid)) +
                     10 + // max length of pos
                     strlen(bcf_f->line->d.id) +
                     strlen(bcf_f->line->d.allele[0]) +
                     strlen(bcf_f->line->d.allele[1]) +
                     4; //tabs
        char *md = (char *) malloc(len * sizeof(char));

        sprintf(md, "%s\t%d\t%s\t%s\t%s",
                     bcf_hdr_id2name(bcf_f->hdr, bcf_f->line->rid),
                     bcf_f->line->pos + 1,
                     bcf_f->line->d.id,
                     bcf_f->line->d.allele[0],
                     bcf_f->line->d.allele[1]); 

        // Write metadata
        md_i += strlen(md);
        md_index[i] = md_i;
        fprintf(md_of, "%s", md);

        // Get gentotypes
        uint32_t num_gts_per_sample = bcf_get_genotypes(bcf_f->hdr,
                                                        bcf_f->line,
                                                        &gt_p,
                                                        &ntmp);
        num_gts_per_sample /= num_inds;
        int32_t *gt_i = gt_p;
        
        // Pack genotypes
        for (j = 0; j < num_inds; ++j) {
            uint32_t gt = 0;
            for (k = 0; k < num_gts_per_sample; ++k) {
                gt += bcf_gt_allele(gt_i[k]);
            }

            packed_ints[int_i] += gt << (30 - 2*two_bit_i);
            two_bit_i += 1;
            if (two_bit_i == 16) {
                two_bit_i = 0;
                int_i += 1;
            }

            sum += gt;
            gt_i += num_gts_per_sample;
        }

        // Get a priority for the variant based on the sum and number of 
        // leading zeros
        p.sum = sum;
        uint32_t prefix_len = 0;
        j = 0;
        while ((j < num_ind_ints) && (packed_ints[j] == 0)){
            prefix_len += 32;
            j += 1;
        }
        if (j < num_ind_ints)
            prefix_len += nlz1(packed_ints[j]);
        
        // Push it into the q
        p.len = prefix_len;
        int *j = (int *) malloc (sizeof(int));
        j[0] = i;
        priq_push(*q, j, p);

        // Write to file
        fwrite(packed_ints, sizeof(uint32_t), num_ind_ints, gt_of);

        memset(packed_ints, 0, num_ind_ints*sizeof(uint32_t));

        free(md);
    }

    free(packed_ints);
    fclose(gt_of);
    fclose(md_of);
}
//}}}

//{{{void sort_gt_md(pri_queue *q,
void sort_gt_md(pri_queue *q,
                uint32_t *md_index,
                uint32_t num_inds,
                uint32_t num_vars,
                char *gt_of_name,
                char *s_gt_of_name,
                char *md_of_name,
                char *bim_out)
{
    FILE *md_of = fopen(md_of_name,"r");
    FILE *md_out = fopen(bim_out,"w");
    FILE *gt_of = fopen(gt_of_name,"rb");
    FILE *s_gt_of = fopen(s_gt_of_name,"wb");

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t *packed_ints = (uint32_t *)
            malloc(num_ind_ints*sizeof(uint32_t));

    priority p;

    // Get variants in order and rewrite a variant-major sorted matrix
    while ( priq_top(*q, &p) != NULL ) {
        int *d = priq_pop(*q, &p);

        uint32_t start = 0;
        if (*d != 0)
            start = md_index[*d - 1];

        uint32_t len = md_index[*d] - start;

        fseek(md_of, start*sizeof(char), SEEK_SET);
        char buf[len+1];
        fread(buf, sizeof(char), len, md_of);
        buf[len] = '\0';

        fseek(gt_of, (*d)*num_ind_ints*sizeof(uint32_t), SEEK_SET);
        fread(packed_ints, sizeof(uint32_t), num_ind_ints, gt_of);
        fwrite(packed_ints, sizeof(uint32_t), num_ind_ints,s_gt_of);

        fprintf(md_out, "%s\n", buf);
    }

    free(packed_ints);

    fclose(md_out);
    fclose(md_of);
    fclose(gt_of);
    fclose(s_gt_of);
}
//}}}


//{{{void rotate_encode_wahbm(uint32_t num_inds,
void rotate_encode_wahbm(uint32_t num_inds,
                         uint32_t num_vars,
                         char *s_gt_of_name,
                         char *r_s_gt_of_name)
{
    uint32_t num_var_ints = 1 + ((num_vars - 1) / 16);
    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t i, j, v;

    uint32_t *I_data = (uint32_t *)calloc(num_var_ints*16,sizeof(uint32_t));
    uint32_t **I = (uint32_t **)malloc(16*sizeof(uint32_t*));
    for (i = 0; i < 16; ++i)
        I[i] = I_data + i*num_var_ints;
    uint32_t I_i = 0, I_int_i = 0, two_bit_i;

    FILE *s_gt_of = fopen(s_gt_of_name,"rb");
    FILE *rs_gt_of = fopen(r_s_gt_of_name,"wb");

    // Write these to values to that this is a well-formed uncompressed 
    // packed int binary file (ubin) file
    //fwrite(&num_vars, sizeof(uint32_t), 1, rs_gt_of);
    //fwrite(&num_inds, sizeof(uint32_t), 1, rs_gt_of);
     
    uint32_t num_inds_to_write = num_inds;
    for (i = 0; i < num_ind_ints; ++i) { // loop over each int col
        for (j = 0; j < num_vars; ++j) { // loop over head row in that col
            // skip to the value at the row/col
            fseek(s_gt_of, 
                  j*num_ind_ints*sizeof(uint32_t) + //row
                  i*sizeof(uint32_t), //col
                  SEEK_SET);

            fread(&v, sizeof(uint32_t), 1, s_gt_of);

            // one int corresponds to a col of 16 two-bit values
            // two_bit_i will move across the cols
            for (two_bit_i = 0; two_bit_i < 16; ++two_bit_i) {
                I[two_bit_i][I_i] += ((v >> (30 - 2*two_bit_i)) & 3) << 
                                     (30 - 2*I_int_i);
            }
            I_int_i += 1;

            if (I_int_i == 16) {
                I_i += 1;
                I_int_i = 0;
            }
        }

        // When we are at the end of the file, and the number of lines 
        // is not a factor of 16, only write out the lines that contain values
        if (num_inds_to_write >= 16) {
            fwrite(I_data,
                   sizeof(uint32_t),
                   num_var_ints*16,
                   rs_gt_of);
            num_inds_to_write -= 16;
        } else {
            fwrite(I_data,
                   sizeof(uint32_t),
                   num_var_ints*num_inds_to_write,
                   rs_gt_of);
        }
        memset(I_data, 0, num_var_ints*16*sizeof(uint32_t));
        I_int_i = 0;
        I_i = 0;
    }

    free(I_data);
    free(I);
    fclose(s_gt_of);
    fclose(rs_gt_of);
}
//}}}

#if 0
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
#endif

//{{{void close_bcf_file(struct bcf_file *bcf_f)
void close_bcf_file(struct bcf_file *bcf_f)
{
    bcf_hdr_destroy(bcf_f->hdr);
    bcf_destroy1(bcf_f->line);
    hts_close(bcf_f->fp);
}
//}}}
