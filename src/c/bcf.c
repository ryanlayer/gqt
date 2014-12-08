#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <assert.h>
#include <inttypes.h>
#include <zlib.h>
#include "genotq.h"
#include "timer.h"


#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

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

//{{{int convert_file_by_name_bcf_to_wahbm_bim(char *in,
int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out,
                                          char *vid_out,
                                          char *tmp_dir)
{
    uint32_t num_inds = num_fields;
    uint32_t num_vars = num_records;

    /*
    char *gt_of_name = ".gt.tmp.packed";
    char *s_gt_of_name = ".s.gt.tmp.packed";
    char *r_s_gt_of_name = ".r.s.gt.tmp.packed";
    char *md_of_name = ".md.tmp.packed";
    char *bim_of_name = ".md.tmp.packed.bim";
    */

    char *gt_of_name,
         *s_gt_of_name,
         *r_s_gt_of_name,
         *md_of_name,
         *bim_of_name;

    asprintf(&gt_of_name, "%s/.gt.tmp.packed", tmp_dir);
    asprintf(&s_gt_of_name, "%s/.s.gt.tmp.packed", tmp_dir);
    asprintf(&r_s_gt_of_name, "%s/.r.s.gt.tmp.packed", tmp_dir);
    asprintf(&md_of_name, "%s/.md.tmp.packed", tmp_dir);
    asprintf(&bim_of_name, "%s/.md.tmp.packed.bim", tmp_dir);


    struct bcf_file bcf_f = init_bcf_file(in);
    pri_queue q = priq_new(0);

    uint64_t *md_index = (uint64_t *) malloc(num_vars * sizeof(uint64_t));

    uint64_t *md_lens = (uint64_t *) malloc(num_vars * sizeof(uint64_t));

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    sort_gt_md(&q,
               md_index,
               md_lens,
               num_inds,
               num_vars,
               gt_of_name,
               s_gt_of_name,
               md_of_name,
               bim_of_name,
               vid_out);

    compress_md(&bcf_f,
                bim_of_name,
                bim_out,
                md_lens,
                num_vars);

    rotate_gt(num_inds,
              num_vars,
              s_gt_of_name,
              r_s_gt_of_name);

    close_bcf_file(&bcf_f);

    int r = convert_file_by_name_ubin_to_wahbm(r_s_gt_of_name, wah_out);

    /*
    remove(gt_of_name);
    remove(s_gt_of_name);
    remove(r_s_gt_of_name);
    remove(md_of_name);
    remove(bim_of_name);
    */


    free(md_index);
    free(md_lens);
    return r;
}
//}}}

//{{{ void push_bcf_gt_md(pri_queue *q,
void push_bcf_gt_md(pri_queue *q,
                    struct bcf_file *bcf_f,
                    uint64_t *md_index,
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

    uint64_t md_i = 0;

    uint32_t i, j, k, sum, int_i, two_bit_i = 0;
    int ntmp = 0;

    int32_t *gt_p = NULL;

    kstring_t md = {0,0,0};

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    fprintf(stderr,"Extracting genotyes and metadata");

    for (i = 0; i < num_vars; ++i) {
        if ((tenth_num_var ==0) || (i % tenth_num_var == 0))
            fprintf(stderr,".");



        sum = 0;
        int_i = 0;
        two_bit_i = 0;

        // Get the next bcf record
        int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);

        // Unpack all of the fields 
        bcf_unpack(bcf_f->line, BCF_UN_ALL);


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

            assert(num_gts_per_sample <= 2);
            assert(num_gts_per_sample > 0);

            if (num_gts_per_sample == 1) {
                if (bcf_gt_is_missing(gt_i[0]))
                    gt = 3;
                else if (bcf_gt_allele(gt_i[k]) == 0)
                    gt = 0;
                else
                    gt = 1;
            } else {
                if (bcf_gt_is_missing(gt_i[0]) || bcf_gt_is_missing(gt_i[1]))
                    gt = 3;
                else if ((bcf_gt_allele(gt_i[0]) == 0 ) &&
                         (bcf_gt_allele(gt_i[1]) == 0 ))
                    gt = 0;
                else if (bcf_gt_allele(gt_i[0]) == bcf_gt_allele(gt_i[1]))
                    gt = 2;
                else if (bcf_gt_allele(gt_i[0]) != bcf_gt_allele(gt_i[1]))
                    gt = 1;

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

        // Get metadata
        bcf_f->line->n_sample = 0;
        vcf_format1(bcf_f->hdr, bcf_f->line, &md);
        md_i += md.l;
        md_index[i] = md_i;
        fprintf(md_of, "%s", md.s);
        md.l = 0;
        memset(packed_ints, 0, num_ind_ints*sizeof(uint32_t));
    }

    fprintf(stderr,"Done\n");

    if (md.s != 0)
        free(md.s);

    free(packed_ints);
    fclose(gt_of);
    fclose(md_of);
}
//}}}

//{{{void sort_gt_md(pri_queue *q,
void sort_gt_md(pri_queue *q,
                uint64_t *md_index,
                uint64_t *md_lens,
                uint32_t num_inds,
                uint32_t num_vars,
                char *gt_of_name,
                char *s_gt_of_name,
                char *md_of_name,
                char *bim_of_name,
                char *vid_out)
{
    // unsorted metadata
    FILE *md_of = fopen(md_of_name,"r");
    // unsorted genotypes
    FILE *gt_of = fopen(gt_of_name,"rb");
    // sorted genotypes
    FILE *md_out = fopen(bim_of_name,"w");
    // sorted variant row #s
    FILE *v_out = fopen(vid_out,"wb");
    // sorted genotypes
    FILE *s_gt_of = fopen(s_gt_of_name,"wb");

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t *packed_ints = (uint32_t *)
            malloc(num_ind_ints*sizeof(uint32_t));

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    uint32_t var_i = 0;
    uint64_t cumul_len = 0;
    fprintf(stderr,"Sorting genotypes and metadata");

    // Get variants in order and rewrite a variant-major sorted matrix
    while ( priq_top(*q, &p) != NULL ) {
        //status
        if ((tenth_num_var == 0) || (var_i % tenth_num_var == 0))
            fprintf(stderr,".");

        // get the data (line num) from the top element
        int *d = priq_pop(*q, &p);

        // get start offset of metadata
        uint64_t start = 0;
        if (*d != 0)
            start = md_index[*d - 1];

        // get the length of meta data
        uint64_t len = md_index[*d] - start;

        // jump to the sport the metadata and read
        fseek(md_of, start*sizeof(char), SEEK_SET);
        char buf[len+1];
        int r = fread(buf, sizeof(char), len, md_of);
        buf[len] = '\0';
        //write metadata to file
        fprintf(md_out, "%s", buf);

        cumul_len += strlen(buf);
        md_lens[var_i] = cumul_len;

        // jump to the sport in the genotypes, read and write
        fseek(gt_of, (*d)*num_ind_ints*sizeof(uint32_t), SEEK_SET);
        r = fread(packed_ints, sizeof(uint32_t), num_ind_ints, gt_of);
        fwrite(packed_ints, sizeof(uint32_t), num_ind_ints,s_gt_of);


        //write out the variant ID
        fwrite(d, sizeof(uint32_t), 1, v_out);

        var_i += 1;
    }

    fprintf(stderr, "Done\n");

    free(packed_ints);

    fclose(md_out);
    fclose(v_out);
    fclose(md_of);
    fclose(gt_of);
    fclose(s_gt_of);
}
//}}}

//{{{ void compress_md(struct bcf_file *bcf_f,
void compress_md(struct bcf_file *bcf_f,
                 char *md_of_name,
                 char *bim_out,
                 uint64_t *md_lens,
                 uint32_t num_vars)
{
    fprintf(stderr, "Compressing metadata.");

    // Get the BCF header
    int h_len;
    char *h_buf = bcf_hdr_fmt_text(bcf_f->hdr, 0, &h_len);

    /* 
     * This header will inclue all of the sample names, which we don't want
     * so we need to scan back to the start of the last line of the header,
     * then clip off the FORMAT and subsequent sample IDs
     */
    //scan back to the start of the header line
    while (h_len && h_buf[h_len-1] != '\n')
        --h_len;
    --h_len;
    while (h_len && h_buf[h_len-1] != '\n') --h_len;
    h_buf[h_len] = '\0';
    --h_len;

    uint64_t c_size = 0;
    uint64_t h_size = strlen(h_buf);
    uint64_t u_size = h_size;
    uint64_t h_i = 0;

    FILE *fp_o = fopen(bim_out, "wb");

     /* 
     * The file is :
     * uncompressed size     ( sizeof(uint64_t))
     * compressed size       ( sizeof(uint64_t))
     * header size           ( sizeof(uint64_t))
     * number of var/records ( bcf_f->num_records*sizeof(uint64_t))
     * md line lengths       ( bcf_f->num_records*sizeof(uint64_t))
     * compressed data 
     */
    fwrite(&u_size, sizeof(uint64_t), 1, fp_o);
    fwrite(&c_size, sizeof(uint64_t), 1, fp_o);
    fwrite(&h_size, sizeof(uint64_t), 1, fp_o);
    uint64_t numv_64 = num_vars;
    fwrite(&numv_64, sizeof(uint64_t), 1, fp_o);
    fwrite(md_lens, sizeof(uint64_t), num_vars, fp_o);

    // in_buf will hold the uncompressed data and out_buf the compressed
    unsigned char *in_buf = (unsigned char *)
            malloc(sizeof(unsigned char) * CHUNK);
    // The output buffer needs to be slightly larger than the in_buff
    unsigned char *out_buf = (unsigned char *)
            malloc(sizeof(unsigned char) * (CHUNK * 2));

    uint32_t in_buf_len = CHUNK;
    uint32_t in_buf_i = 0;
    uint32_t out_buf_len = CHUNK * 2;

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    int ret = deflateInit(&strm, 6);
    if (ret != Z_OK) {
        fprintf(stderr, "error: Cannot init stream\n");
        exit(1);
    }

    // have is used to store the size of the compressed data
    uint64_t have;

    //{{{ Compress header
    while (h_i < h_size) {
        /* 
         * Move either the size of the buffer, or the remaining amout of the
         * header  to the buffer
         */
        uint32_t mv_len = MIN(h_size - h_i, in_buf_len - in_buf_i);
        memcpy(in_buf + in_buf_i, h_buf + h_i, mv_len);

        // Move the head pointer of both the buffer and f1
        in_buf_i += mv_len;
        h_i += mv_len;

        /* 
         * When the buffer is full, compress it to out_buf and write it to a
         * file
         */
        if (in_buf_i == in_buf_len) {
            strm.avail_in = in_buf_len;
            strm.next_in = in_buf;

            strm.avail_out = out_buf_len;
            strm.next_out = out_buf;

            ret = deflate(&strm, Z_FULL_FLUSH);

            if (ret == Z_BUF_ERROR) {
                fprintf(stderr,
                        "No progress is possible; either avail_in or "
                        "avail_out was zero %u\t%u.\n", 
                        strm.avail_in, strm.avail_out);
                exit(1);
            } else if (ret == Z_MEM_ERROR) {
                fprintf(stderr, "Insufficient memory.\n");
                exit(1);
            } else if (ret == Z_STREAM_ERROR) {
                fprintf(stderr,
                        "The state (as represented in stream) is inconsistent,"
                        " or stream was NULL.");
                exit(1);
            }

            // The amount compressed is the amount of the buffer used
            have = out_buf_len - strm.avail_out;

            // Track the size of the compressed data
            c_size += have;
            if (fwrite(out_buf, 1, have, fp_o) != have) {
                fprintf(stderr, "Error writing compressed value 0.\n");
                exit(1);
            }
            //fwrite(in_buf, 1, CHUNK, fp_o);
            in_buf_i = 0;
        }
    }
    //}}}

    // Read and compress the sorted meta data fields into the buffer
    FILE *fp = NULL;

    struct stat md_stat;
    stat(md_of_name, &md_stat);
    size_t md_size = md_stat.st_size;
    uint64_t md_size_tenth = md_size/10;
    uint64_t md_size_status = md_size_tenth;

    if ((fp = fopen(md_of_name, "r")) == NULL) {
        fprintf(stderr,
                "error: unable to open file %s.\n",
                md_of_name);
        exit(1);
    }

    /* 
     * It is likely that there is still data on the buffer to be compressed.
     * to_read is the amount left, start by reading and compressing that
     * then move on the whole chunks
     */
    uint32_t read = 0, to_read = in_buf_len - in_buf_i;


    while ((read = fread(in_buf + in_buf_i, sizeof(char), to_read, fp)) ==
            to_read ) {
        u_size += read;

        strm.next_in = in_buf;
        strm.avail_in = in_buf_len;
        strm.next_out = out_buf;
        strm.avail_out = out_buf_len;

        ret = deflate(&strm, Z_FULL_FLUSH);

        if (ret == Z_BUF_ERROR) {
            fprintf(stderr,
                    "No progress is possible; either avail_in or "
                    "avail_out was zero %u\t%u.\n", 
                    strm.avail_in, strm.avail_out);
            exit(1);
        } else if (ret == Z_MEM_ERROR) {
            fprintf(stderr, "Insufficient memory.\n");
            exit(1);
        } else if (ret == Z_STREAM_ERROR) {
            fprintf(stderr,
                    "The state (as represented in stream) is inconsistent,"
                    " or stream was NULL.");
            exit(1);
        }

        have = out_buf_len - strm.avail_out;

        // Track the size of the compressed data
        c_size += have;
        if (fwrite(out_buf, 1, have, fp_o) != have) {
            fprintf(stderr, "Error writing compressed value 1.\n");
            exit(1);
        }

        in_buf_i = 0;
        to_read = in_buf_len;

        if ( u_size > md_size_status) {
            md_size_status += md_size_tenth;
            fprintf(stderr, ".");
        }
    }

    // It is likely that there is still data on the buffer to be compressed.

    if (read > 0) {
        u_size += read;

        //fwrite(in_buf, 1, in_buf_i + read, fp_o);
        strm.next_in = in_buf;
        strm.avail_in = in_buf_i + read;
        strm.next_out = out_buf;
        strm.avail_out = out_buf_len;
        ret = deflate(&strm, Z_FULL_FLUSH);

        if (ret == Z_BUF_ERROR) {
            fprintf(stderr,
                    "No progress is possible; either avail_in or "
                    "avail_out was zero %u\t%u.\n", 
                    strm.avail_in, strm.avail_out);
            exit(1);
        } else if (ret == Z_MEM_ERROR) {
            fprintf(stderr, "Insufficient memory.\n");
            exit(1);
        } else if (ret == Z_STREAM_ERROR) {
            fprintf(stderr,
                    "The state (as represented in stream) is inconsistent,"
                    " or stream was NULL.");
            exit(1);
        }

        have = out_buf_len - strm.avail_out;

        // Track the size of the compressed data
        c_size += have;
        if (fwrite(out_buf, 1, have, fp_o) != have) {
            fprintf(stderr, "Error writing compressed value 1.\n");
            exit(1);
        }
    }


    // update the header values
    fseek(fp_o, 0, SEEK_SET);
    fwrite(&u_size, sizeof(uint64_t), 1, fp_o);
    fwrite(&c_size, sizeof(uint64_t), 1, fp_o);

    fclose(fp_o);
    fclose(fp);
    free(in_buf);
    free(out_buf);
    free(h_buf);

    fprintf(stderr, "Done\n");
}
//}}}

//{{{void rotate_gt(uint32_t num_inds,
void rotate_gt(uint32_t num_inds,
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
    fwrite(&num_vars, sizeof(uint32_t), 1, rs_gt_of);
    fwrite(&num_inds, sizeof(uint32_t), 1, rs_gt_of);
     
    uint32_t tenth_num_ind_ints = num_ind_ints / 10;
    fprintf(stderr, "Rotating genotypes");


    uint32_t num_inds_to_write = num_inds;
    for (i = 0; i < num_ind_ints; ++i) { // loop over each int col
        if ((tenth_num_ind_ints == 0) || (i % tenth_num_ind_ints == 0))
            fprintf(stderr,".");

        for (j = 0; j < num_vars; ++j) { // loop over head row in that col
            // skip to the value at the row/col
            fseek(s_gt_of, 
                  j*num_ind_ints*sizeof(uint32_t) + //row
                  i*sizeof(uint32_t), //col
                  SEEK_SET);

            int r = fread(&v, sizeof(uint32_t), 1, s_gt_of);

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
    fprintf(stderr,"Done\n");

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
#endif

//{{{void close_bcf_file(struct bcf_file *bcf_f)
void close_bcf_file(struct bcf_file *bcf_f)
{
    bcf_hdr_destroy(bcf_f->hdr);
    bcf_destroy1(bcf_f->line);
    hts_close(bcf_f->fp);
}
//}}}
