#define _GNU_SOURCE
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <assert.h>
#include <inttypes.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sysexits.h>

#include "bcf.h"
#include "bm.h"
#include "wah.h"
#include "variant_metadata.h"
#include "vid.h"
#include "off.h"
#include "bim.h"
#include "ubin.h"
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
    
    bcf_f.file_name = file_name;

    bcf_f.fp.bcf = hts_open(file_name,"rb");
    if ( !bcf_f.fp.bcf ) 
        err(EX_DATAERR, "Could not read file: %s", file_name);

    if (bcf_f.fp.bcf->format.compression != bgzf)
        errx(EX_DATAERR, "Not a BGZF file: %s\n", file_name);

    bcf_f.hdr = bcf_hdr_read(bcf_f.fp.bcf);

    if ( !bcf_f.hdr )
        err(EX_DATAERR, "Could not read the header: %s", file_name);

    htsFormat type = *hts_get_format(bcf_f.fp.bcf);
    if (type.format == bcf) {
        bcf_f.is_bcf = 1;
        bcf_f.offset = bgzf_tell(bcf_f.fp.bcf->fp.bgzf);
    } else {
        bcf_f.is_bcf = 0;

        bcf_f.str.m = 0;
        bcf_f.str.l = 0;
        bcf_f.str.s = 0;

        hts_close(bcf_f.fp.bcf);
        bcf_f.fp.vcf = bgzf_open(file_name, "r");
        
        if (!bcf_f.fp.vcf)
            err(EX_DATAERR, "Could not read file: %s", file_name);
        if ( !(bcf_f.fp.vcf->is_compressed) ) {
            bgzf_close(bcf_f.fp.vcf);
            err(EX_DATAERR, "Not a compresed file: %s", file_name);
        }

        // move past the header
        int ret;
        uint64_t last_offset = 0;
        while (bgzf_getline(bcf_f.fp.vcf, '\n', &(bcf_f.str)) >= 0) {
            if (bcf_f.str.s[0] == '#')
                last_offset = bgzf_tell(bcf_f.fp.vcf);
            else
                break;
        }
        if (bgzf_seek(bcf_f.fp.vcf, last_offset, SEEK_SET) != 0)
            err(EX_DATAERR, "Error moving past header: %s", file_name);
        bcf_f.offset = last_offset;
    }

    bcf_f.line = bcf_init1();
    bcf_f.num_records = bcf_hdr_nsamples(bcf_f.hdr);
    bcf_f.gt = NULL;

    return bcf_f;
}
//}}}

//{{{int get_bcf_line(struct bcf_file *bcf_f)
int get_bcf_line(struct bcf_file *bcf_f)
{
    int r;

    if (bcf_f->is_bcf) {
        r = bcf_read(bcf_f->fp.bcf, bcf_f->hdr, bcf_f->line);
        bcf_f->offset = bgzf_tell(bcf_f->fp.bcf->fp.bgzf);
    } else {
        r = bgzf_getline(bcf_f->fp.vcf, '\n', &(bcf_f->str));
        if (r < 0)
            return r;
        else 
            r = 0;
        bcf_f->offset = bgzf_tell(bcf_f->fp.vcf);
        vcf_parse(&(bcf_f->str), bcf_f->hdr, bcf_f->line);
    }

    return r;
}
//}}}

//{{{ int goto_bcf_line(struct bcf_file *bcf_f,
int goto_bcf_line(struct bcf_file *bcf_f,
                   struct off_file *off_f,
                   uint32_t line_no)
{
    int r;

    if (line_no >= off_f->gqt_header->num_variants)
        errx(EX_SOFTWARE,
             "Seeking to line beyond file boundary for '%s'.",
             off_f->file_name); 


    if (bcf_f->is_bcf) 
        r = bgzf_seek(bcf_f->fp.bcf->fp.bgzf,
                      off_f->offsets[line_no],
                      SEEK_SET);
    else
        r = bgzf_seek(bcf_f->fp.vcf,
                      off_f->offsets[line_no],
                      SEEK_SET);

    return r;
}
//}}}

//{{{int convert_file_by_name_bcf_to_wahbm_offset(char *in,
int convert_file_by_name_bcf_to_wahbm_offset(char *in,
                                             uint32_t num_fields,
                                             uint32_t num_records,
                                             char *wah_out,
                                             char *offset_out,
                                             char *vid_out,
                                             char *tmp_dir,
                                             char *full_cmd)
{
    uint32_t num_inds = num_fields;
    uint32_t num_vars = num_records;

    char *gt_of_name,
         *gt_s_of_name,
         *gt_s_r_of_name;

    int r = asprintf(&gt_of_name, "%s/.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_of_name, "%s/.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_r_of_name, "%s/.r.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");


    struct bcf_file bcf_f = init_bcf_file(in);
    pri_queue q = priq_new(0);

    push_bcf_gt_offset(&q,
                       &bcf_f,
                       num_inds,
                       num_vars,
                       gt_of_name,
                       offset_out,
                       full_cmd);
    sort_gt(&q,
            num_inds,
            num_vars,
            gt_of_name,
            gt_s_of_name,
            vid_out,
            full_cmd);

    rotate_gt(num_inds,
              num_vars,
              gt_s_of_name,
              gt_s_r_of_name);

    close_bcf_file(&bcf_f);

    r = convert_file_by_name_ubin_to_wahbm(gt_s_r_of_name,
                                          wah_out,
                                          full_cmd);

    remove(gt_of_name);
    remove(gt_s_of_name);
    remove(gt_s_r_of_name);

    return r;
}
//}}}

//{{{int convert_file_by_name_bcf_to_wahbm_bim(char *in,
int convert_file_by_name_bcf_to_wahbm_bim(char *in,
                                          uint32_t num_fields,
                                          uint32_t num_records,
                                          char *wah_out,
                                          char *bim_out,
                                          char *vid_out,
                                          char *tmp_dir,
                                          char *full_cmd)
{
    uint32_t num_inds = num_fields;
    uint32_t num_vars = num_records;

    char *gt_of_name,
         *gt_s_of_name,
         *gt_s_r_of_name,
         *md_of_name,
         *md_s_of_name;

    int r = asprintf(&gt_of_name, "%s/.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_of_name, "%s/.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_r_of_name, "%s/.r.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&md_of_name, "%s/.md.tmp", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&md_s_of_name, "%s/.s.md.tmp", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");


    struct bcf_file bcf_f = init_bcf_file(in);
    pri_queue q = priq_new(0);

    uint64_t *md_index = (uint64_t *) malloc(num_vars * sizeof(uint64_t));
    if (!md_index )
        err(EX_OSERR, "malloc error");
    
    uint64_t *md_s_index = (uint64_t *) malloc(num_vars * sizeof(uint64_t));
    if (!md_s_index )
        err(EX_OSERR, "malloc error");

    push_bcf_gt_md(&q,
                   &bcf_f,
                   md_index,
                   num_inds,
                   num_vars,
                   gt_of_name,
                   md_of_name);

    sort_gt(&q,
            num_inds,
            num_vars,
            gt_of_name,
            gt_s_of_name,
            vid_out,
            full_cmd);

    compress_md(&bcf_f,
                md_of_name,
                bim_out,
                md_index,
                num_vars,
                num_inds,
                full_cmd);

    rotate_gt(num_inds,
              num_vars,
              gt_s_of_name,
              gt_s_r_of_name);

    close_bcf_file(&bcf_f);

    r = convert_file_by_name_ubin_to_wahbm(gt_s_r_of_name,
                                               wah_out,
                                               full_cmd);

    remove(gt_of_name);
    remove(gt_s_of_name);
    remove(gt_s_r_of_name);
    remove(md_of_name);
    remove(md_s_of_name);

    //free(md_index);
    free(md_s_index);
    return r;
}
//}}}

//{{{int convert_file_by_name_bcf_to_wahbm_metadata_offset(char *in,
int convert_file_by_name_bcf_to_wahbm_metadata_offset(char *in,
                                                      uint32_t num_fields,
                                                      uint32_t num_records,
                                                      char *wah_out,
                                                      char *bim_out,
                                                      char *offset_out,
                                                      char *vid_out,
                                                      char *tmp_dir,
                                                      char *full_cmd)
{
    uint32_t num_inds = num_fields;
    uint32_t num_vars = num_records;

    char *gt_of_name,
         *gt_s_of_name,
         *gt_s_r_of_name,
         *md_of_name,
         *md_s_of_name;

    int r = asprintf(&gt_of_name, "%s/.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_of_name, "%s/.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&gt_s_r_of_name, "%s/.r.s.gt.tmp.packed", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&md_of_name, "%s/.md.tmp", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");
    r = asprintf(&md_s_of_name, "%s/.s.md.tmp", tmp_dir);
    if (r == -1) err(EX_OSERR, "asprintf error");


    struct bcf_file bcf_f = init_bcf_file(in);
    pri_queue q = priq_new(0);

    uint64_t *md_index = (uint64_t *) malloc(num_vars * sizeof(uint64_t));
    if (!md_index )
        err(EX_OSERR, "malloc error");
    
    uint64_t *md_s_index = (uint64_t *) malloc(num_vars * sizeof(uint64_t));
    if (!md_s_index )
        err(EX_OSERR, "malloc error");


    push_bcf_gt_md_offset(&q,
                          &bcf_f,
                          md_index,
                          num_inds,
                          num_vars,
                          gt_of_name,
                          md_of_name,
                          offset_out,
                          full_cmd);
    sort_gt(&q,
            num_inds,
            num_vars,
            gt_of_name,
            gt_s_of_name,
            vid_out,
            full_cmd);

    compress_md(&bcf_f,
                md_of_name,
                bim_out,
                md_index,
                num_vars,
                num_inds,
                full_cmd);

    rotate_gt(num_inds,
              num_vars,
              gt_s_of_name,
              gt_s_r_of_name);

    close_bcf_file(&bcf_f);

    r = convert_file_by_name_ubin_to_wahbm(gt_s_r_of_name,
                                               wah_out,
                                               full_cmd);

    remove(gt_of_name);
    remove(gt_s_of_name);
    remove(gt_s_r_of_name);
    remove(md_of_name);
    remove(md_s_of_name);

    //free(md_index);
    free(md_s_index);
    return r;
}
//}}}

//{{{ void push_bcf_gt_offset(pri_queue *q,
void push_bcf_gt_offset(pri_queue *q,
                       struct bcf_file *bcf_f,
                       uint32_t num_inds,
                       uint32_t num_vars,
                       char *gt_of_name,
                       char *offset_of_name,
                       char *full_cmd)
{
    FILE *gt_of = fopen(gt_of_name,"wb");
    if (!gt_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", gt_of_name);

    /*
    FILE *off_of = fopen(offset_of_name,"wb");
    if (!off_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", offset_of_name);
    */

    struct off_file *off_of = new_off_file(offset_of_name,
                                           full_cmd,
                                           num_vars,
                                           num_inds);

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    // Allocate this array once, and reuse it for every variant
    uint32_t *packed_ints = (uint32_t *) calloc(num_ind_ints,
                                                sizeof(uint32_t));

    uint64_t md_i = 0;

    uint32_t i, j, k, sum, int_i, two_bit_i = 0;
    int ntmp = 0;

    int32_t *gt_p = NULL;

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    fprintf(stderr,"Extracting genotypes and offsets");

    uint64_t last_offset = 0;

    for (i = 0; i < num_vars; ++i) {
        if ((tenth_num_var ==0) || (i % tenth_num_var == 0))
            fprintf(stderr,".");

        sum = 0;
        int_i = 0;
        two_bit_i = 0;

        last_offset = bcf_f->offset;
        // Get the next bcf record
        int r = get_bcf_line(bcf_f);

        if (r == -1) 
            err(EX_NOINPUT, "Error reading file \"%s\"", bcf_f->file_name);
        
        // Unpack all of the fields 
        bcf_unpack(bcf_f->line, BCF_UN_ALL);


        // Get gentotypes
        uint32_t num_gts_per_sample = bcf_get_genotypes(bcf_f->hdr,
                                                        bcf_f->line,
                                                        &gt_p,
                                                        &ntmp);
        num_gts_per_sample /= num_inds;
        int32_t *gt_i = gt_p;

        
        if (num_gts_per_sample != 2) {
            fprintf(stderr, "num_gts_per_sample:%u\t%u:%u\n",
                            num_gts_per_sample,
                            i,
                            num_vars);
        }
        // Pack genotypes
        for (j = 0; j < num_inds; ++j) {
            uint32_t gt = 0;

            assert(num_gts_per_sample <= 2);
            assert(num_gts_per_sample > 0);

            if ( (num_gts_per_sample == 1) || 
                 (gt_i[1] == bcf_int32_vector_end) ){
                if (bcf_gt_is_missing(gt_i[0]))
                    gt = 3;
                else if (bcf_gt_allele(gt_i[0]) == 0)
                    gt = 0;
                else
                    gt = 1;
            } else {
                if (bcf_gt_is_missing(gt_i[0]) && bcf_gt_is_missing(gt_i[1]))
                    gt = 3;
                else if ((bcf_gt_allele(gt_i[0]) == 0 ) &&
                         (bcf_gt_allele(gt_i[1]) == 0 ))
                    gt = 0;
                else if ((bcf_gt_allele(gt_i[0]) != (bcf_gt_allele(gt_i[1]))))
                    gt = 1;
                else 
                    gt = 2;
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
        if (!j)
            err(EX_OSERR, "malloc error");
        j[0] = i;
        priq_push(*q, j, p);

        // Write to file
        if (fwrite(packed_ints,
                   sizeof(uint32_t),
                   num_ind_ints,
                   gt_of) != num_ind_ints)
            err(EX_IOERR, "Error writing to \"%s\"", gt_of_name); 

        // Write offset the of the start of the current line
        add_to_off_file(off_of, last_offset);

        memset(packed_ints, 0, num_ind_ints*sizeof(uint32_t));
    }

    fprintf(stderr,"Done\n");

    free(packed_ints);
    fclose(gt_of);
    destroy_off_file(off_of);
}
//}}}

//{{{ void push_bcf_gt_md(pri_queue *q,
void push_bcf_gt_md_offset(pri_queue *q,
                           struct bcf_file *bcf_f,
                           uint64_t *md_index,
                           uint32_t num_inds,
                           uint32_t num_vars,
                           char *gt_of_name,
                           char *md_of_name,
                           char *offset_of_name,
                           char *full_cmd)
{

    FILE *gt_of = fopen(gt_of_name,"wb");
    if (!gt_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", gt_of_name);


    FILE *md_of = fopen(md_of_name,"w");
    if (!md_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", md_of_name);

    struct off_file *off_of = new_off_file(offset_of_name,
                                           full_cmd,
                                           num_vars,
                                           num_inds);

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    // Allocate this array once, and reuse it for every variant
    uint32_t *packed_ints = (uint32_t *) calloc(num_ind_ints,
                                                sizeof(uint32_t));

    uint64_t md_i = 0;

    uint32_t i, j, k, sum, int_i, two_bit_i = 0;
    int ntmp = 0;

    int32_t *gt_p = NULL;

    kstring_t md = {0,0,0};

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    fprintf(stderr,"Extracting genotypes and metadata");

    uint64_t last_offset = 0;

    for (i = 0; i < num_vars; ++i) {
        if ((tenth_num_var ==0) || (i % tenth_num_var == 0))
            fprintf(stderr,".");

        sum = 0;
        int_i = 0;
        two_bit_i = 0;

        // Get the next bcf record
        //int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);
        //int r = bcf_read(bcf_f->fp.bcf, bcf_f->hdr, bcf_f->line);

        last_offset = bcf_f->offset;

        int r = get_bcf_line(bcf_f);

        if (r == -1) 
            err(EX_NOINPUT, "Error reading file \"%s\"", bcf_f->file_name);
        
        // Unpack all of the fields 
        bcf_unpack(bcf_f->line, BCF_UN_ALL);


        // Get gentotypes
        uint32_t num_gts_per_sample = bcf_get_genotypes(bcf_f->hdr,
                                                        bcf_f->line,
                                                        &gt_p,
                                                        &ntmp);
        num_gts_per_sample /= num_inds;
        int32_t *gt_i = gt_p;

        
        if (num_gts_per_sample != 2) {
            fprintf(stderr, "num_gts_per_sample:%u\t%u:%u\n",
                            num_gts_per_sample,
                            i,
                            num_vars);
        }
        // Pack genotypes
        for (j = 0; j < num_inds; ++j) {
            uint32_t gt = 0;

            assert(num_gts_per_sample <= 2);
            assert(num_gts_per_sample > 0);

            if ( (num_gts_per_sample == 1) || 
                 (gt_i[1] == bcf_int32_vector_end) ){
                if (bcf_gt_is_missing(gt_i[0]))
                    gt = 3;
                else if (bcf_gt_allele(gt_i[0]) == 0)
                    gt = 0;
                else
                    gt = 1;
            } else {
                if (bcf_gt_is_missing(gt_i[0]) && bcf_gt_is_missing(gt_i[1]))
                    gt = 3;
                else if ((bcf_gt_allele(gt_i[0]) == 0 ) &&
                         (bcf_gt_allele(gt_i[1]) == 0 ))
                    gt = 0;
                else if ((bcf_gt_allele(gt_i[0]) != (bcf_gt_allele(gt_i[1]))))
                    gt = 1;
                else 
                    gt = 2;
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
        if (!j)
            err(EX_OSERR, "malloc error");
        j[0] = i;
        priq_push(*q, j, p);

        // Write to file
        if (fwrite(packed_ints,
                   sizeof(uint32_t),
                   num_ind_ints,
                   gt_of) != num_ind_ints)
            err(EX_IOERR, "Error writing to \"%s\"", gt_of_name); 

        // Get metadata
        bcf_f->line->n_sample = 0;
        vcf_format1(bcf_f->hdr, bcf_f->line, &md);
        md_i += md.l;
        md_index[i] = md_i;
        fprintf(md_of, "%s", md.s);
        md.l = 0;

        add_to_off_file(off_of, last_offset);

        memset(packed_ints, 0, num_ind_ints*sizeof(uint32_t));
    }

    fprintf(stderr,"Done\n");

    if (md.s != 0)
        free(md.s);

    free(packed_ints);
    fclose(gt_of);
    fclose(md_of);
    destroy_off_file(off_of);
}
//}}}

//{{{ void push_bcf_gt_md(pri_queue *q,
/*
 * Read a BCF file and populated:
 * INPUT:
 *   bcf_f: a bcf_file struct that contains the file handle and metadata for
 *          the target BCF file 
 * OUTPUT:
 *   q: a prioriy queue where the priority is the alt allele count and
 *      the number of leaded zeros, and the value is the index of that 
 *      variant both in md_of_name (where md_index maps index to the file
 *      offset) and in gt_of_name (which is fixed width and that index can be
 *      multiplied by width to get fill offset)
 *   md_of_name: a file containg the variant metadata
 *   md_index: an array of 64-bit ints that contains that maps the index of a
 *             vairiant to its offset in md_of_name
 *   gt_of_name: a variant-major file containing arrays of pakced genotypes
 *               (0,1,2,3) this is a fixed-width file and each element in the
 *               array is a 32-bit int
 */
void push_bcf_gt_md(pri_queue *q,
                    struct bcf_file *bcf_f,
                    uint64_t *md_index,
                    uint32_t num_inds,
                    uint32_t num_vars,
                    char *gt_of_name,
                    char *md_of_name)
{

    FILE *gt_of = fopen(gt_of_name,"wb");
    if (!gt_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", gt_of_name);


    FILE *md_of = fopen(md_of_name,"w");
    if (!md_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", md_of_name);

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    // Allocate this array once, and reuse it for every variant
    uint32_t *packed_ints = (uint32_t *) calloc(num_ind_ints,
                                                sizeof(uint32_t));

    uint64_t md_i = 0;

    uint32_t i, j, k, sum, int_i, two_bit_i = 0;
    int ntmp = 0;

    int32_t *gt_p = NULL;

    kstring_t md = {0,0,0};

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    fprintf(stderr,"Extracting genotypes and metadata");

    for (i = 0; i < num_vars; ++i) {
        if ((tenth_num_var ==0) || (i % tenth_num_var == 0))
            fprintf(stderr,".");

        sum = 0;
        int_i = 0;
        two_bit_i = 0;

        // Get the next bcf record
        //int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);
        //int r = bcf_read(bcf_f->fp.bcf, bcf_f->hdr, bcf_f->line);
        int r = get_bcf_line(bcf_f);

        if (r == -1) 
            err(EX_NOINPUT, "Error reading file \"%s\"", bcf_f->file_name);
        
        // Unpack all of the fields 
        bcf_unpack(bcf_f->line, BCF_UN_ALL);


        // Get gentotypes
        uint32_t num_gts_per_sample = bcf_get_genotypes(bcf_f->hdr,
                                                        bcf_f->line,
                                                        &gt_p,
                                                        &ntmp);
        num_gts_per_sample /= num_inds;
        int32_t *gt_i = gt_p;

        
        if (num_gts_per_sample != 2) {
            fprintf(stderr, "num_gts_per_sample:%u\t%u:%u\n",
                            num_gts_per_sample,
                            i,
                            num_vars);
        }
        // Pack genotypes
        for (j = 0; j < num_inds; ++j) {
            uint32_t gt = 0;

            assert(num_gts_per_sample <= 2);
            assert(num_gts_per_sample > 0);

            if ( (num_gts_per_sample == 1) || 
                 (gt_i[1] == bcf_int32_vector_end) ){
                if (bcf_gt_is_missing(gt_i[0]))
                    gt = 3;
                else if (bcf_gt_allele(gt_i[0]) == 0)
                    gt = 0;
                else
                    gt = 1;
            } else {
                if (bcf_gt_is_missing(gt_i[0]) && bcf_gt_is_missing(gt_i[1]))
                    gt = 3;
                else if ((bcf_gt_allele(gt_i[0]) == 0 ) &&
                         (bcf_gt_allele(gt_i[1]) == 0 ))
                    gt = 0;
                else if ((bcf_gt_allele(gt_i[0]) != (bcf_gt_allele(gt_i[1]))))
                    gt = 1;
                else 
                    gt = 2;
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
        if (!j)
            err(EX_OSERR, "malloc error");
        j[0] = i;
        priq_push(*q, j, p);

        // Write to file
        if (fwrite(packed_ints,
                   sizeof(uint32_t),
                   num_ind_ints,
                   gt_of) != num_ind_ints)
            err(EX_IOERR, "Error writing to \"%s\"", gt_of_name); 

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
/*
 * Use the priority queue q to create new variant metadata and packed int files
 * that are sorted.  Also create a md_s_index, which stores the end of each
 * variant metadata entry as part of the bim file.
 *
 * Each element is poped from the q.  The value of that element (d) is used to
 * find the variants in gt_of_name (d*width) and the metadata inf md_of_name
 * (md_index[d-1] or 0 if d is zero).  The target variant line is written to
 * gt_s_of_name and the target metadata is writein to md_s_of_name.  The length
 * of that metadata line is added with the previous and appended to md_s_index.
 * The index (d) is appended to vid_out.
 *
 * INPUT:
 *   q: a prioriy queue where the priority is the alt allele count and
 *      the number of leaded zeros, and the value is the index of that 
 *      variant both in md_of_name (where md_index maps index to the file
 *      offset) and in gt_of_name (which is fixed width and that index can be
 *      multiplied by width to get fill offset)
 *   md_of_name: a file containg the variant metadata
 *   md_index: an array of 64-bit ints that contains that maps the index of a
 *             vairiant to its offset in md_of_name
 *   gt_of_name: a variant-major file containing arrays of pakced genotypes
 *               (0,1,2,3) this is a fixed-width file and each element in the
 *               array is a 32-bit int
 *
 * OUTPUT:
 *   md_s_index:  an array of 64-bit ints that store the end of each metadata
 *                entry in the sorted md_s_of_name file
 *   gt_s_of_name: a variant-major file containing sorted arrays of pakced
 *                 genotypes (0,1,2,3) this is a fixed-width file and each
 *                 element in the array is a 32-bit int
 *   md_s_of_name: a file containing the sorted variant metadata
 *   vid_out: a file containg the index of the sorted variant in the original
 *            BCF file
 */
void sort_gt(pri_queue *q,
             uint32_t num_inds,
             uint32_t num_vars,
             char *gt_of_name,
             char *gt_s_of_name,
             char *vid_out,
             char *full_cmd)
{
    // unsorted genotypes
    FILE *gt_of = fopen(gt_of_name,"rb");
    if (!gt_of)
        err(EX_NOINPUT, "Cannot read file\"%s\"", gt_of_name);

    struct vid_file *v_out = new_vid_file(vid_out,
                                          full_cmd,
                                          num_vars,
                                          num_inds); 
 
    // sorted genotypes
    FILE *s_gt_of = fopen(gt_s_of_name,"wb");
    if (!s_gt_of)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", gt_s_of_name);

    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t *packed_ints = (uint32_t *)
            malloc(num_ind_ints*sizeof(uint32_t));
    if (!packed_ints )
        err(EX_OSERR, "malloc error");

    priority p;

    uint32_t tenth_num_var = num_vars / 10;
    uint32_t var_i = 0;
    uint64_t cumul_len = 0;
    fprintf(stderr,"Sorting genotypes");

    // Get variants in order and rewrite a variant-major sorted matrix
    while ( priq_top(*q, &p) != NULL ) {
        //status
        if ((tenth_num_var == 0) || (var_i % tenth_num_var == 0))
            fprintf(stderr,".");

        // get the data (line num) from the top element
        int *d = priq_pop(*q, &p);

        uint64_t start = 0;
        int r;

        // jump to the sport in the genotypes, read and write
        start = num_ind_ints*sizeof(uint32_t);
        start = (*d)*start;
        fseek(gt_of, start, SEEK_SET);

        size_t fr = fread(packed_ints, sizeof(uint32_t), num_ind_ints, gt_of);
        check_file_read(gt_of_name, gt_of, num_ind_ints, fr);

        if (fwrite(packed_ints,
                   sizeof(uint32_t),
                   num_ind_ints,
                   s_gt_of) != num_ind_ints)
            err(EX_IOERR, "Error writing to \"%s\"", gt_s_of_name); 

        write_vid(v_out, *d);

        var_i += 1;
    }

    fprintf(stderr, "Done\n");

    free(packed_ints);

    destroy_vid_file(v_out);
    fclose(gt_of);
    fclose(s_gt_of);
}
//}}}

//{{{ void compress_md(struct bcf_file *bcf_f,
/*
 * Create the bim file that contains some header information then the
 * compressed metadata. md line lengths are from md_s_index
 *
 * The file is :
 * uncompressed size     ( sizeof(uint64_t))
 * compressed size       ( sizeof(uint64_t))
 * header size           ( sizeof(uint64_t))
 * number of var/records ( bcf_f->num_records*sizeof(uint64_t))
 * md line lengths       ( bcf_f->num_records*sizeof(uint64_t))
 * compressed data 
 *
 * INPUT
 *   bcf_f: a bcf_file struct that contains the file handle and metadata for
 *          the target BCF file 
 *   md_s_of_name: a file containing the sorted variant metadata
 *   md_s_index:  an array of 64-bit ints that store the end of each metadata
 *                entry in the sorted md_s_of_name file
 *   num_var: number of variants
 * OUTPUT
 *   bim_out: the file containing some metadata about the variant metadata and
 *            the compressed metadata
 */
void compress_md(struct bcf_file *bcf_f,
                 char *md_s_of_name,
                 char *bim_out,
                 uint64_t *md_s_index,
                 uint32_t num_vars,
                 uint32_t num_inds,
                 char *full_cmd)
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

#if 1
    struct bim_file *bim_f = new_bim_file(bim_out,
                                          full_cmd,
                                          num_vars,
                                          num_inds,
                                          u_size,
                                          c_size,
                                          h_size,
                                          md_s_index);
#endif

#if 0
    FILE *fp_o = fopen(bim_out, "wb");
    if (!fp_o)
        err(EX_CANTCREAT, "Cannot create file \"%s\"", bim_out);

     /* 
     * The file is :
     * uncompressed size     ( sizeof(uint64_t))
     * compressed size       ( sizeof(uint64_t))
     * header size           ( sizeof(uint64_t))
     * number of var/records ( bcf_f->num_records*sizeof(uint64_t))
     * md line lengths       ( bcf_f->num_records*sizeof(uint64_t))
     * compressed data 
     */
    if (fwrite(&u_size, sizeof(uint64_t), 1, fp_o) != 1)
            err(EX_IOERR, "Error writing to \"%s\"", bim_out); 

    if (fwrite(&c_size, sizeof(uint64_t), 1, fp_o) != 1)
            err(EX_IOERR, "Error writing to \"%s\"", bim_out); 

    if (fwrite(&h_size, sizeof(uint64_t), 1, fp_o) != 1)
            err(EX_IOERR, "Error writing to \"%s\"", bim_out); 

    uint64_t numv_64 = num_vars;

    if (fwrite(&numv_64, sizeof(uint64_t), 1, fp_o) != 1)
            err(EX_IOERR, "Error writing to \"%s\"", bim_out); 

    if (fwrite(md_s_index, sizeof(uint64_t), num_vars, fp_o) != num_vars)
            err(EX_IOERR, "Error writing to \"%s\"", bim_out); 
#endif

    // in_buf will hold the uncompressed data and out_buf the compressed
    unsigned char *in_buf = (unsigned char *)
            malloc(sizeof(unsigned char) * CHUNK);
    if (!in_buf )
        err(EX_OSERR, "malloc error");

    // The output buffer needs to be slightly larger than the in_buff
    unsigned char *out_buf = (unsigned char *)
            malloc(sizeof(unsigned char) * (CHUNK * 2));
    if (!out_buf)
        err(EX_OSERR, "malloc error");

    uint32_t in_buf_len = CHUNK;
    uint32_t in_buf_i = 0;
    uint32_t out_buf_len = CHUNK * 2;

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    int ret = deflateInit(&strm, 6);
    if (ret != Z_OK) {
        errx(EX_CANTCREAT,"Zlib error: Cannot init stream\n");
    }

    // have is used to store the size of the compressed data
    uint64_t have;

    seek_bim_to_data(bim_f);

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
                errx(EX_SOFTWARE,
                        "Zlib error: No progress is possible; either avail_in or "
                        "avail_out was zero %u\t%u.", 
                        strm.avail_in, strm.avail_out);
            } else if (ret == Z_MEM_ERROR) {
                errx(EX_SOFTWARE, "Zlib error: Insufficient memory.");
            } else if (ret == Z_STREAM_ERROR) {
                errx(EX_SOFTWARE,
                        "Zlib error: The state (as represented in stream) is "
                        "inconsistent or stream was NULL.");
            }

            // The amount compressed is the amount of the buffer used
            have = out_buf_len - strm.avail_out;

            // Track the size of the compressed data
            c_size += have;
            //if (fwrite(out_buf, 1, have, fp_o) != have) 
            assert(bim_f->type == BIM_LOCAL);
            if (fwrite(out_buf, 1, have, bim_f->file.local) != have) 
                err(EX_IOERR, "Error writing compressed value 0.");

            //fwrite(in_buf, 1, CHUNK, fp_o);
            in_buf_i = 0;
        }
    }
    //}}}

    // Read and compress the sorted meta data fields into the buffer
    FILE *fp = NULL;

    struct stat md_stat;
    stat(md_s_of_name, &md_stat);
    size_t md_size = md_stat.st_size;
    uint64_t md_size_tenth = md_size/10;
    uint64_t md_size_status = md_size_tenth;

    fp = fopen(md_s_of_name, "r");
    if (!fp)
        err(EX_NOINPUT, "Cannot read file \"%s\"", md_s_of_name);

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
            errx(EX_SOFTWARE,
                    "Zlib error: No progress is possible; either avail_in or "
                    "avail_out was zero %u\t%u.", 
                    strm.avail_in, strm.avail_out);
        } else if (ret == Z_MEM_ERROR) {
            errx(EX_SOFTWARE, "Zlib error: Insufficient memory.");
        } else if (ret == Z_STREAM_ERROR) {
            errx(EX_SOFTWARE,
                 "Zlib error: The state (as represented in stream) is "
                 "inconsistent or stream was NULL.");
        }

        have = out_buf_len - strm.avail_out;

        // Track the size of the compressed data
        c_size += have;
        //if (fwrite(out_buf, 1, have, fp_o) != have) 
        if (bim_f->type != BIM_LOCAL)
            errx(1, "Cannot write to remote BIM file.");

        if (fwrite(out_buf, 1, have, bim_f->file.local) != have) 
            err(EX_IOERR, "Error writing compressed value 1.");

        in_buf_i = 0;
        to_read = in_buf_len;

        if ( u_size > md_size_status) {
            md_size_status += md_size_tenth;
            fprintf(stderr, ".");
        }
    }

    if (ferror(fp))
        err(EX_IOERR, "Error reading file \"%s\"", md_s_of_name);

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
            errx(EX_SOFTWARE,
                 "Zlib error: No progress is possible; either avail_in or "
                 "avail_out was zero %u\t%u.", 
                 strm.avail_in, strm.avail_out);
        } else if (ret == Z_MEM_ERROR) {
            errx(EX_SOFTWARE, "Zlib error: Insufficient memory.");
        } else if (ret == Z_STREAM_ERROR) {
            errx(EX_SOFTWARE,
                 "Zlib error: The state (as represented in stream) is "
                 "inconsistent or stream was NULL.");
        }

        have = out_buf_len - strm.avail_out;

        // Track the size of the compressed data
        c_size += have;
        //if (fwrite(out_buf, 1, have, fp_o) != have) 
        if (bim_f->type != BIM_LOCAL)
            errx(1,"Cannot write to remote BIM file.");
        if (fwrite(out_buf, 1, have, bim_f->file.local) != have) 
            err(EX_IOERR, "Error writing compressed value 1.");
    }


#if 0
    // update the header values
    fseek(fp_o, 0, SEEK_SET);
    if (fwrite(&u_size, sizeof(uint64_t), 1, fp_o) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", bim_out); 

    if (fwrite(&c_size, sizeof(uint64_t), 1, fp_o) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", bim_out); 
    fclose(fp_o);
#endif

    update_bim_file_header(u_size, c_size, h_size, bim_f);

    destroy_bim_file(bim_f);
    
    fclose(fp);
    free(in_buf);
    free(out_buf);
    free(h_buf);

    fprintf(stderr, "Done\n");
}
//}}}

//{{{void rotate_gt(uint32_t num_inds,
/* Take a variant major packed integer genotype file and rotate it to be a
 * individual major file.
 *
 * In the rotated file, the first variant (row) in the variant major file will
 * become the first column in the individual major file.  The rotation is done
 * blocked manner.  That is,  1 int (16 genotypes) will be read from a line in
 * the variant-major file.  Those genotypes will be appended to a matrix of 16
 * individual-major ints.  We then jump to the same column in the next variant,
 * read another int (16 genotypes), and then transpose those into the
 * individual-major matrix.  After 16 variants have been read, that matrix is
 * written to the file, and the process continues, with the next group of 16
 * variants.  Once the first int (16 genotypes) have been read and transposed,
 * we return to the next int (next 16 genotypes) in the first 16 variants and
 * proceed similarly.
 *
 * INPUT
 *   num_inds: number of individual
 *   num_var: number of variants
 *   gt_s_of_name: a variant-major file containing sorted arrays of pakced
 *                 genotypes (0,1,2,3) this is a fixed-width file and each
 *                 element in the array is a 32-bit int
 * OUTPUT
 *   gt_s_r_of_name: an individual-major file containing sorted arrays of
 *                   pakced genotypes (0,1,2,3) this is a fixed-width file and
 *                   each element in the array is a 32-bit int
 */
void rotate_gt(uint32_t num_inds,
               uint32_t num_vars,
               char *gt_s_of_name,
               char *gt_s_r_of_name)
{
    uint32_t num_var_ints = 1 + ((num_vars - 1) / 16);
    uint32_t num_ind_ints = 1 + ((num_inds - 1) / 16);

    uint32_t i, j, v;

    uint32_t *I_data = (uint32_t *)calloc(num_var_ints*16,sizeof(uint32_t));
    if (!I_data )
        err(EX_OSERR, "calloc error");

    uint32_t **I = (uint32_t **)malloc(16*sizeof(uint32_t*));
    if (!I)
        err(EX_OSERR, "malloc error");

    for (i = 0; i < 16; ++i)
        I[i] = I_data + i*num_var_ints;
    uint32_t I_i = 0, I_int_i = 0, two_bit_i;

    FILE *s_gt_of = fopen(gt_s_of_name,"rb");
    if (!s_gt_of)
        err(EX_NOINPUT, "Cannot read file\"%s\"", gt_s_of_name);

    FILE *rs_gt_of = fopen(gt_s_r_of_name,"wb");
    if (!rs_gt_of)
        err(EX_NOINPUT, "Cannot read file\"%s\"", gt_s_r_of_name);

    // Write these to values to a well-formed uncompressed packed int binary
    // file (ubin) file
    if (fwrite(&num_vars, sizeof(uint32_t), 1, rs_gt_of) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", gt_s_r_of_name); 
    if (fwrite(&num_inds, sizeof(uint32_t), 1, rs_gt_of) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", gt_s_r_of_name); 
     
    uint32_t tenth_num_ind_ints = num_ind_ints / 10;
    fprintf(stderr, "Rotating genotypes");


    uint32_t num_inds_to_write = num_inds;
    for (i = 0; i < num_ind_ints; ++i) { // loop over each int col
        if ((tenth_num_ind_ints == 0) || (i % tenth_num_ind_ints == 0))
            fprintf(stderr,".");

        for (j = 0; j < num_vars; ++j) { // loop over head row in that col
            // skip to the value at the row/col
            uint64_t row = j;
            row *=  num_ind_ints;
            row *=  sizeof(uint32_t);
            uint64_t col = i;
            col *= sizeof(uint32_t);
            fseek(s_gt_of, row + col, SEEK_SET);

            size_t fr = fread(&v, sizeof(uint32_t), 1, s_gt_of);
            check_file_read(gt_s_of_name, s_gt_of, 1, fr);

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
            if (fwrite(I_data,
                   sizeof(uint32_t),
                   num_var_ints*16,
                   rs_gt_of) != num_var_ints*16)
                err(EX_IOERR, "Error writing to \"%s\"", gt_s_r_of_name); 

            num_inds_to_write -= 16;
        } else {
            if (fwrite(I_data,
                   sizeof(uint32_t),
                   num_var_ints*num_inds_to_write,
                   rs_gt_of) != num_var_ints*num_inds_to_write)
                err(EX_IOERR, "Error writing to \"%s\"", gt_s_r_of_name); 
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

//{{{void close_bcf_file(struct bcf_file *bcf_f)
void close_bcf_file(struct bcf_file *bcf_f)
{
    bcf_hdr_destroy(bcf_f->hdr);
    bcf_destroy1(bcf_f->line);
    hts_close(bcf_f->fp);
}
//}}}

//{{{uint32_t get_variant_metadata(struct bcf_file *bcf_f,
uint32_t get_variant_metadata(struct bcf_file *bcf_f,
                              uint32_t num_vars,
                              char *field_name,
                              void **values,
                              void *missing_value,
                              int *type)
{

    float **f_val = (float **) values;
    //int **i_val = (int **) values;

    float *f_missing = (float *)missing_value;
    //int *i_missing = (int *)missing_value;

    int id = bcf_hdr_id2int(bcf_f->hdr,BCF_DT_ID,field_name);
    if ( !bcf_hdr_idinfo_exists(bcf_f->hdr,BCF_HL_INFO,id) ) {
        fprintf(stderr,
                "The INFO tag \"%s\" is not defined in the header\n",
                field_name);
        return 1;
    }

    *type = bcf_hdr_id2type(bcf_f->hdr,BCF_HL_INFO,id);

    if (*type == BCF_HT_FLAG)
        //*i_val = (int *) malloc(num_vars * sizeof(int));
        *f_val = (float *) malloc(num_vars * sizeof(float));
    else if (*type == BCF_HT_INT)
        //*i_val = (int *) malloc(num_vars * sizeof(int));
        *f_val = (float *) malloc(num_vars * sizeof(float));
    else if (*type == BCF_HT_REAL)
        *f_val = (float *) malloc(num_vars * sizeof(float));
    else {
        fprintf(stderr,
                "The type for INFO tag \"%s\" is not supported\n",
                field_name);
        return 1;
    }

    uint32_t i;

    for (i = 0; i < num_vars; ++i) {
        // Get the next bcf record
        int r = bcf_read(bcf_f->fp, bcf_f->hdr, bcf_f->line);
        r = bcf_hdr_set_samples(bcf_f->hdr, NULL, 0);
        
        bcf_unpack(bcf_f->line, BCF_UN_INFO);
        bcf_info_t *info = bcf_get_info(bcf_f->hdr, bcf_f->line, field_name);

        if ( (*type != BCF_HT_FLAG) && (info->len != 1) ) {
            fprintf(stderr,
                    "Non-scalar value found in field \"%s\" in variant at"
                    "%d:%d\n", 
                    field_name,
                    bcf_f->line->rid,
                    bcf_f->line->pos);
            return 1;
        }

        if (info == NULL) {
            if (*type == BCF_HT_FLAG)
                (*f_val)[i] = 0;
            else if (*type == BCF_HT_INT)
                (*f_val)[i] = *f_missing;
            else if (*type == BCF_HT_REAL)
                (*f_val)[i] = *f_missing;
        } else {
            if (*type == BCF_HT_FLAG)
                (*f_val)[i] = 1;
            else if (*type == BCF_HT_INT) {
                (*f_val)[i] = (float)info->v1.i;
            } else if (*type == BCF_HT_REAL)
                (*f_val)[i] = info->v1.f;
        }
    }

    return 0;
}
//}}}

//{{{ int get_variant_metadata_type(struct bcf_file *bcf_f,
int get_variant_metadata_type(struct bcf_file *bcf_f,
                              char *field_name)
{
    int id = bcf_hdr_id2int(bcf_f->hdr,BCF_DT_ID,field_name);
    return bcf_hdr_id2type(bcf_f->hdr,BCF_HL_INFO,id);
}
//}}}

//{{{ uint32_t index_variant_metadata(char *bcf_file_name,
uint32_t index_variant_metadata(char *bcf_file_name,
                                char *vid_file_name,
                                char *db_file_name,
                                char *variant_metadata_index_file_name,
                                uint32_t num_variants,
                                char *field_name,
                                uint32_t num_to_test,
                                uint32_t *num_bins,
                                void **bin_range_lo,
                                void **bin_range_hi,
                                int *less_than_bin,
                                int *greater_than_bin)
{
    if (num_to_test > num_variants)
        num_to_test = num_variants;

    struct bcf_file bcf_f = init_bcf_file(bcf_file_name);

    int type = get_variant_metadata_type(&bcf_f, field_name);

    void *missing;
    int int_missing = 0;
    float float_missing = 0.0;

    if (type == BCF_HT_FLAG) {
        *less_than_bin = 0;
        *greater_than_bin = 0;
    } else {
        if (*less_than_bin == -1) {
            *less_than_bin = 1;
        }

        if (*greater_than_bin == -1) {
            *greater_than_bin = 1;
        }
    }

    if (type == BCF_HT_REAL) {
        missing = &float_missing;
    } else if ((type == BCF_HT_INT) || (type == BCF_HT_FLAG)) {
        missing = &float_missing;
    }

    void *variant_metadata;

    int r = get_variant_metadata(&bcf_f,
                                 num_variants,
                                 field_name,
                                 &variant_metadata,
                                 missing,
                                 &type);

    close_bcf_file(&bcf_f);

    if ((type == BCF_HT_INT) || 
        (type == BCF_HT_REAL) ||
        (type == BCF_HT_FLAG)) {

        float *float_variant_metadata = (float *)variant_metadata;
        float **float_bin_range_lo = (float **)bin_range_lo,
              **float_bin_range_hi = (float **)bin_range_hi;

         uint32_t actual_num_bins =
                float_equal_freq_binning(float_variant_metadata,
                                         num_to_test,
                                         *num_bins,
                                         float_bin_range_lo,
                                         float_bin_range_hi);
     
         *num_bins = actual_num_bins + *less_than_bin + *greater_than_bin;

#if 0
        // open VID file
        FILE *vid_f = fopen(vid_file_name, "rb");
        if (vid_f == NULL) {
            fprintf(stderr, "Could not read VID file: %s\n", vid_file_name);
            return 1;
        }
        uint32_t *vids = (uint32_t *) malloc(num_variants*sizeof(uint32_t));
        r = fread(vids, sizeof(uint32_t), num_variants, vid_f);
        fclose(vid_f);
#endif
        
        struct vid_file *vf = open_vid_file(vid_file_name);
        load_vid_data(vf);

        /* Since the variants are in allele freq order, we need to copy
         * the resulting value to an array that is back in the original
         * order
         */
        float *mapped_float_variant_metadata = 
                (float *)calloc(num_variants, sizeof(float));
        uint32_t i;
        for ( i = 0; i < num_variants; ++i) {
            mapped_float_variant_metadata[i] = 
                float_variant_metadata[vf->vids[i]];
        }

        destroy_vid_file(vf);
        free(float_variant_metadata);
        float_variant_metadata = mapped_float_variant_metadata;



        uint32_t **ws = NULL;
        uint32_t *w_lens = NULL;

        uint32_t bit_array_lens = float_to_wah_bitmap(float_variant_metadata,
                                                      num_variants,
                                                      *float_bin_range_lo,
                                                      *float_bin_range_hi,
                                                      actual_num_bins,
                                                      *less_than_bin,
                                                      *greater_than_bin,
                                                      *num_bins,
                                                      &ws,
                                                      &w_lens); 

        FILE *f = fopen(variant_metadata_index_file_name, "a");

        if (!f) {
            fprintf(stderr,
                    "Unable to open %s\n",
                    variant_metadata_index_file_name);
            return 1;
        }

        long start_offset = ftell(f);

        fwrite(w_lens, sizeof(uint32_t), *num_bins, f);

        for ( i = 0; i < *num_bins; ++i) {
            fwrite(ws[i], sizeof(uint32_t), w_lens[i], f);
            free(ws[i]);
        }

        fclose(f);

        int rowid = register_variant_metadata_index(bcf_file_name,
                                                    db_file_name,
                                                    field_name,
                                                    start_offset,
                                                    type,
                                                    *num_bins);
        if (rowid == -1) {
            exit(1);
        }

        int rows_added = add_variant_metadata_float_bins(db_file_name,
                                                         rowid,
                                                         *float_bin_range_lo,
                                                         *float_bin_range_hi,
                                                         actual_num_bins,
                                                         *less_than_bin,
                                                         *greater_than_bin);
    }

    return type;
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
<<<<<<< HEAD
=======
>>>>>>> no_bim
=======
>>>>>>> f3afedbc6485986b0a86fa9cab7cb65c4ae3f7fe

    if (bcf_f->is_bcf) {
        hts_close(bcf_f->fp.bcf);
    } else {
        if (bcf_f->str.s)
            free(bcf_f->str.s);
        bgzf_close(bcf_f->fp.vcf);
    }
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
    if (!md )
        err(EX_OSERR, "malloc error");

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
