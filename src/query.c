#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <inttypes.h>
#include "genotq.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

int query_help();

void print_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        uint32_t *vids,
                        struct gqt_query *q,
                        uint32_t **counts,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *bim);
int query_cmp(uint32_t value,
              int op_condition,
              int condition_value);


void get_bcf_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        struct gqt_query *q,
                        char **id_query_list,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *vid_file_name,
                        char *src_bcf_file_name,
                        int bcf_output);

int compare_uint32_t (const void *a, const void *b);

//{{{ int compare_uint32_t (const void *a, const void *b)
int compare_uint32_t (const void *a, const void *b)
{
    return ( *(uint32_t*)a - *(uint32_t*)b );
}
//}}}

//{{{int query_cmp(uint32_t value,
int query_cmp(uint32_t value,
              int op_condition,
              int condition_value)
{
    switch(op_condition) {
        case p_equal:
            return value == condition_value;
        case p_not_equal:
            return value != condition_value;
        case p_less_than:
            return value < condition_value;
        case p_greater_than:
            return value > condition_value;
        case p_less_than_equal:
            return value <= condition_value;
        case p_greater_than_equal:
            return value >= condition_value;
        default:
            return -1;
    }
}
//}}}

//{{{ int query(int argc, char **argv)
int query(int argc, char **argv)
{
    if (argc < 2) return query_help();

    int c;
    char *wahbm_file_name=NULL,
         *id_query=NULL,
         *gt_query=NULL,
         *db_file_name=NULL,
         *bim_file_name=NULL,
         *src_bcf_file_name=NULL,
         *vid_file_name=NULL;
    int i_is_set = 0,
        id_q_count = 0,
        gt_q_count = 0,
        d_is_set = 0,
        c_is_set = 0,
        v_is_set = 0,
        s_is_set = 0,
        b_is_set = 0,
        bcf_output = 0;

    char *id_query_list[100];
    char *gt_query_list[100];

    //{{{ parse cmd line opts
    while ((c = getopt (argc, argv, "chi:p:g:d:b:v:s:B")) != -1) {
        switch (c) {
        case 'c':
            c_is_set = 1;
            break;
        case 'i':
            i_is_set = 1;
            wahbm_file_name = optarg;
            break;
        case 'p':
            id_query_list[id_q_count] = optarg;
            id_q_count += 1;
            break;
        case 'g':
            gt_query_list[gt_q_count] = optarg;
            gt_q_count += 1;
            break;
        case 'd':
            d_is_set = 1;
            db_file_name = optarg;
            break;
        case 'b':
            b_is_set = 1;
            bim_file_name = optarg;
            break;
        case 'v':
            v_is_set = 1;
            vid_file_name = optarg;
            break;
        case 's':
            s_is_set = 1;
            src_bcf_file_name = optarg;
            break;
        case 'B':
            bcf_output = 1;
            break;
        case 'h':
            return query_help();
        case '?':
            if ( (optopt == 'i') ||
                    (optopt == 'p') ||
                    (optopt == 'g') ||
                    (optopt == 'd') ||
                    (optopt == 'b') )
                fprintf (stderr, "Option -%c requires an argument.\n",
                         optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        default:
            return query_help();
        }
    }

    // Try to auto-detect file names based on GQT
    if ( (i_is_set == 1) && (b_is_set == 0)) {

        int auto_bim_file_name_size = asprintf(&bim_file_name,
                                               "%s",
                                               wahbm_file_name);
        strcpy(bim_file_name + strlen(bim_file_name) - 3, "bim");

        if ( access( bim_file_name, F_OK) != -1 ) {
            b_is_set = 1;
        } else {
            printf("Auto detect failure: BIM file %s not found\n",
                   bim_file_name);
            return query_help();
        }
    }

    if ( (i_is_set == 1) && (v_is_set == 0)) {

        int auto_vid_file_name_size = asprintf(&vid_file_name,
                                               "%s",
                                               wahbm_file_name);
        strcpy(vid_file_name + strlen(vid_file_name) - 3, "vid");

        if ( access( vid_file_name, F_OK) != -1 ) {
            v_is_set = 1;
        } else {
            printf("Auto detect failure: VID file %s not found\n",
                   vid_file_name);
            return query_help();
        }
    }


#if 0
    // Try to auto-detect file names based on GQT
    if ( (v_is_set == 0) && (i_is_set == 1)) {

        int auto_vid_file_name_size = asprintf(&vid_file_name,
                                               "%s",
                                               wahbm_file_name);

        strcpy(vid_file_name + strlen(vid_file_name) - 3, "vid");

        if ( access( vid_file_name, F_OK) != -1 ) {
            v_is_set = 1;
        } else {
            printf("Auto detect failure: VID file %s not found\n",
                    vid_file_name);
            return query_help();
        }

        int auto_wahbm_file_name_size = asprintf(&wahbm_file_name,
                                                 "%s.gqt",
                                                 src_bcf_file_name);

        if ( access( wahbm_file_name, F_OK) != -1 ) {
            i_is_set = 1;
        } else {
            printf("Auto detect failure: WAH file %s not found\n",
                   wahbm_file_name);
            return query_help();
        }
    }
#endif

#if 0
    if (((v_is_set == 1) && (s_is_set == 0)) ||
        ((v_is_set == 0) && (s_is_set == 1)) ) {
        printf("Either BOTH VID and source BCF must be set, or neither.\n");
        return query_help();
    }

    if (((v_is_set == 1) && (s_is_set == 1)) && (b_is_set == 1)) {
        printf("Set EITHER VID/soruce BCF or bim, not BOTH.\n");
        return query_help();
    } 
 
    if (((v_is_set == 0) && (s_is_set == 0)) && (b_is_set == 0)) {
        printf("Must set EITHER VID/soruce BCF or bim.\n");
        return query_help();
    } 
#endif

    if (i_is_set == 0) {
        printf("GQT file is not set\n");
        return query_help();
    }

    if (v_is_set == 0) {
        printf("VID file is not set\n");
        return query_help();
    }

    if (b_is_set == 0) {
        printf("BIM file is not set\n");
        return query_help();
    }

    if (d_is_set == 0) {
        printf("PED database file is not set\n");
        return query_help();
    }

    if (gt_q_count != id_q_count) {
        printf("Mismatched number of individual and genotype query strings\n");
        return query_help();
    }
    //}}}

    struct gqt_query q[100];
    uint32_t *gt_mask[100];
    uint32_t *counts[100];
    uint32_t *mapped_counts[100];
    uint32_t id_lens[100];

    int r, i, j, k;

    for (i = 0; i < gt_q_count; ++i) {
        if (parse_q(gt_query_list[i], &(q[i]))) {
            fprintf(stderr, "in the %dth genotype query.\n", i+1);
            return 1;
        }
    }

    // open WAH/GQT file
    struct wah_file wf = init_wahbm_file(wahbm_file_name);

    // open VID file
    FILE *vid_f = fopen(vid_file_name, "rb");
    if (vid_f == NULL) {
        fprintf(stderr, "Could not read VIDE file: %s\n", vid_file_name);
        return 1;
    }
    uint32_t *vids = (uint32_t *) malloc(wf.num_fields*sizeof(uint32_t));
    r = fread(vids, sizeof(uint32_t), wf.num_fields, vid_f);
    fclose(vid_f);

    uint32_t num_ints = (wf.num_fields + 32 - 1)/ 32;
    uint32_t len_ints;

    for (i = 0; i < gt_q_count; ++i) {
        uint32_t len_count_R;
        uint32_t *R;
        /* 
         * Submit the population query to the PED database and get back both
         * the list of of ids in R and the length of R in id_lens[i]
         */
        id_lens[i] = resolve_ind_query(&R,
                                      id_query_list[i],
                                      db_file_name);

        // Enforce that the offsets of the relevant samples is 
        // within the number of samples in the GQT index.
        if (id_lens[i] > wf.num_records) {
            fprintf(stderr, 
                    "ERROR: there are more samples in the PED database (%d) "
                    "that match this condition \nthan there are in the GQT "
                    "index (%d).  Perhaps your PED file is a superset of "
                    "the\nsamples in your VCF/BCF file?\n", 
                    id_lens[i], 
                    wf.num_records);
            return 1;
        }

        uint32_t low_v, high_v;

        /*
         * q holds the parameters of each query, first determin the range of 
         * bitmaps to pull
         */
        if ( q[i].variant_op == p_maf ) {
            low_v = 1;
            high_v = 3;
        } else {
            if ( q[i].genotype_condition[0] == 1)
                low_v = 0;
            else if ( q[i].genotype_condition[1] == 1)
                low_v = 1;
            else if ( q[i].genotype_condition[2] == 1)
                low_v = 2;
            else if ( q[i].genotype_condition[3] == 1)
                low_v = 3;

            if ( q[i].genotype_condition[3] == 1)
                high_v = 4;
            else if ( q[i].genotype_condition[2] == 1)
                high_v = 3;
            else if ( q[i].genotype_condition[1] == 1)
                high_v = 2;
            else if ( q[i].genotype_condition[0] == 1)
                high_v = 1;
        }

        /*
         * The set of variants that are printed is stored in a mask for each
         * query, then those masks are combine to a final mask.  Each mask is a
         * 32-bit packed int, where each bit correspons to one variant.  How
         * those bits are set depends on the filter the user specifices.
         *
         * If they simply ask for a count or perecent, then there is not filter
         * and the mask is set to all 1s.
         *
         * If count is followed by a condition, then the count/pct is compared
         * to that condition and the bits are set for those that meet the
         * condition.
         *
         * If no funtion is used then we simply run the wahbm range query and
         * convert the wah results to packed ints for the mask
         *
         */

        /* User asks for a count, percent, or maf */
        if ( ( q[i].variant_op == p_count ) || 
             ( q[i].variant_op == p_pct ) ||
             ( q[i].variant_op == p_maf ) ) {

            if (q[i].variant_op == p_maf) {
#ifdef __AVX2__
            len_count_R = avx_sum_range_records_in_place_wahbm(wf,
                                                               R,
                                                               id_lens[i],
                                                               low_v,
                                                               high_v,
                                                               &(counts[i]));
#else
            len_count_R = sum_range_records_in_place_wahbm(wf,
                                                           R,
                                                           id_lens[i],
                                                           low_v,
                                                           high_v,
                                                           &(counts[i]));
                                                           
#endif
            } else {
#ifdef __AVX2__
                len_count_R = 
                    avx_count_range_records_in_place_wahbm(wf,
                                                           R,
                                                           id_lens[i],
                                                           low_v,
                                                           high_v,
                                                           &(counts[i]));
#else
                len_count_R = 
                    count_range_records_in_place_wahbm(wf,
                                                       R,
                                                       id_lens[i],
                                                       low_v,
                                                       high_v,
                                                       &(counts[i]));
#endif
            }

            /* Since the variants are in allele freq order, we need to copy
             * the resulting value to an array that is back in the original
             * order
             */
            mapped_counts[i] = (uint32_t *)calloc(len_count_R,
                                                  sizeof(uint32_t));
            for ( j = 0; j < len_count_R; ++j)
                mapped_counts[i][vids[j]] = counts[i][j];

            gt_mask[i] = (uint32_t *) malloc(num_ints * sizeof(uint32_t));

            /* User specifies a condition */
            if ( q[i].op_condition != -1) { 

                /* Since we only find counts, when the user asks for a
                 * perecent, just convert that back to the count that meets the
                 * percent condition
                 */
                float condition_value = q[i].condition_value;
                if (q[i].variant_op == p_pct) 
                    condition_value *= id_lens[i];
                else if (q[i].variant_op == p_maf)
                    condition_value *= id_lens[i]*2;


                /* Test to see if each count meets the condition */
                uint32_t v = 0, int_i = 0, bit_i = 0;
                for ( j = 0; j < len_count_R; ++j) {
                    if ( query_cmp(counts[i][j],
                                   q[i].op_condition,
                                   condition_value) ) {
                        v |= 1 << (31 - bit_i);
                    }

                    bit_i += 1;
                    if ( bit_i == 32 ) {
                        gt_mask[i][int_i] = v;
                        int_i += 1;
                        bit_i = 0;
                        v = 0;
                    }
                }
            
                if ( bit_i > 0)
                    gt_mask[i][int_i] = v;
            } else {
                // if no op is set then let everything pass
                for (j = 0; j < num_ints; ++j)
                    gt_mask[i][j] = -1; // set all the bits to 1
            }
        /* User only gives genotype filters, no funtion/condition */
        } else {
            uint32_t *gt_R;
            uint32_t len_wf_R = range_records_in_place_wahbm(wf,
                                                                 R,
                                                                 id_lens[i],
                                                                 low_v,
                                                                 high_v,
                                                                 &gt_R);
            len_ints = wah_to_ints(gt_R,len_wf_R,&(gt_mask[i]));
            free(gt_R);
        }
        free(R);
    }

    uint32_t *final_mask = (uint32_t *) calloc(num_ints,sizeof(uint32_t));

    // combine all of the masks to see what we need to print
    for (i = 0; i < num_ints; ++i) {
        final_mask[i] = ~0;
        for (j = 0; j < gt_q_count; ++j)
            final_mask[i] &= gt_mask[j][i];
    }

    if (c_is_set == 1) {
        uint32_t masked_vid_count = 0;

        for (i = 0; i < num_ints; ++i) 
            masked_vid_count += popcount(final_mask[i]);

        if (masked_vid_count <= wf.num_fields)
            printf("%u\n", masked_vid_count);
        else
            printf("%u\n", wf.num_fields);

    } else if ((v_is_set == 1) && (s_is_set == 1)) {
        get_bcf_query_result(final_mask,
                             num_ints, 
                             q,
                             id_query_list,
                             id_lens,
                             gt_q_count,
                             wf.num_fields,
                             vid_file_name,
                             src_bcf_file_name,
                             bcf_output);
    } else if (b_is_set == 1){

        uint32_t *mapped_mask = (uint32_t *) calloc(num_ints,sizeof(uint32_t));

        uint32_t v,p,leading_zeros, hit;
        for (i = 0; i < num_ints; ++i) {
            if (final_mask[i] != 0) {
                v = final_mask[i];
                p = popcount(v);
                for (j = 0; j < p; ++j) {
                    leading_zeros = __builtin_clz(v);

                    if (i*32 + leading_zeros + 1 > wf.num_fields)
                        break;

                    hit = vids[leading_zeros + i*32];

                    mapped_mask[hit/32] |= 1 << (31-hit%32);
                    v &= ~(1 << (32 - leading_zeros - 1));
                }
            }
            if (i*32 + leading_zeros + 1 > wf.num_fields)
                break;
        }

        print_query_result(mapped_mask,
                           num_ints,
                           vids,
                           q,
                           mapped_counts,
                           id_lens,
                           gt_q_count,
                           wf.num_fields,
                           bim_file_name);
    }

    for (j = 0; j < gt_q_count; ++j) {
        free(gt_mask[j]);
        if ( (q[j].variant_op == p_count) || 
             (q[j].variant_op == p_pct) ||
             (q[j].variant_op == p_maf) )
            free(counts[j]);
    }

    fclose(wf.file);
    return 0;
}
//}}}

//{{{ void get_bcf_query_result(uint32_t *mask,
void get_bcf_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        struct gqt_query *q,
                        char **id_query_list,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *vid_file_name,
                        char *src_bcf_file_name,
                        int bcf_output)
{

    /* The VID file contains the line numbers of the variants after they have
     * been sorted.  To reach back into the BCF file to print the metadata
     * associated with the variants marked in the mask, we need to create a
     * sorted list of line numbers we want.  So first we intersect the VID file
     * and the mask, then sort it.
     */
    FILE *vid_f = fopen(vid_file_name, "rb");
    uint32_t *vids = (uint32_t *) malloc(num_fields*sizeof(uint32_t));
    int r = fread(vids, sizeof(uint32_t), num_fields, vid_f);
    fclose(vid_f);

    uint32_t i, j, masked_vid_count = 0;

    for (i = 0; i < mask_len; ++i)
        masked_vid_count += popcount(mask[i]);

    uint32_t *masked_vids = (uint32_t *)
            malloc(masked_vid_count*sizeof(uint32_t));
    uint32_t masked_vid_i = 0;

    for (i = 0; i < mask_len; ++i) {
        uint32_t bytes = mask[i];
	if (bytes == 0)
            continue; /* skip a bunch of ops if you can */
        for (j = 0; j < 32; j++) {
            if (bytes & (1 << (31 - j))) {
                masked_vids[masked_vid_i] = vids[i*32 + j];
                masked_vid_i+=1;
            }
        }
        if (masked_vid_i == masked_vid_count)
            break;
    }

    free(vids);

    qsort(masked_vids, masked_vid_count, sizeof(uint32_t), compare_uint32_t);

    htsFile *fp    = hts_open(src_bcf_file_name,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *line    = bcf_init1();
    //bcf_hdr_set_samples(hdr, print_name_csv, 0);

    htsFile *out;
    if (!bcf_output)
        out = hts_open("-", "w");
    else
        out = hts_open("-", "wb");

    r = bcf_hdr_write(out, hdr);

    uint32_t bcf_line_i = 0;
    masked_vid_i = 0;
    while ( bcf_read(fp, hdr, line) != -1) {
        if (masked_vids[masked_vid_i] == bcf_line_i) {
            r = bcf_unpack(line, BCF_UN_ALL);
            r = bcf_write1(out, hdr, line);
            masked_vid_i+=1;
        }
        if (masked_vid_i == masked_vid_count)
            break;
        bcf_line_i += 1;
    }
    hts_close(out);
    hts_close(fp);
}
//}}}
 
//{{{ void print_query_result(uint32_t *mask,
void print_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        uint32_t *vids,
                        struct gqt_query *q,
                        uint32_t **counts,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *bim)
{
    uint32_t i,j,k,line_idx,bytes, bit_i = 0;

    struct quick_file_info qfile;
    struct output_buffer outbuf;
    char pct[50];


    init_out_buf(&outbuf, NULL);
    quick_file_init(bim, &qfile);


    append_out_buf(&outbuf,
                   qfile.main_buf,
                   qfile.header_len);

    char *info_s;
    for (k=0; k < num_qs; k++) {
        if ( q[k].variant_op == p_count ) {
            asprintf(&info_s, "##INFO=<ID=GTQ_%u,Number=1,Type=Integer,"
                              "Description=\"GQT count result from query "
                              "%u\">\n",
                              k, k);
            append_out_buf(&outbuf, info_s, strlen(info_s));
        } else if ( q[k].variant_op == p_pct ) {
            asprintf(&info_s, "##INFO=<ID=GTQ_%u,Number=1,Type=Float,"
                              "Description=\"GQT percent result from query "
                              "%u\">\n",
                              k, k);
            append_out_buf(&outbuf, info_s, strlen(info_s));
        } else if ( q[k].variant_op == p_maf ) {
            asprintf(&info_s, "##INFO=<ID=GTQ_%u,Number=1,Type=Float,"
                              "Description=\"GQT maf result from query "
                              "%u\">\n",
                              k, k);
            append_out_buf(&outbuf, info_s, strlen(info_s));
        }

    }

    char last_header_line[]="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    append_out_buf(&outbuf, last_header_line, strlen(last_header_line));



    for (i=0; i < mask_len; ++i) {
        bytes = mask[i];
	if (bytes == 0)
            continue; /* skip a bunch of ops if you can */
        for (j=0; j < 32; j++) {
            if (bytes & 1 << (31 - j)) {
	        line_idx = i*32+j;
	        append_out_buf(&outbuf,
                               qfile.lines[line_idx],
                               qfile.line_lens[line_idx]-1);
                for (k=0; k < num_qs; k++) {
                    if ( q[k].variant_op == p_count ) {
                        asprintf(&info_s, ";GTQ_%u=%u", k,
                                counts[k][line_idx]);
                                //counts[k][vids[line_idx]]);
                        append_out_buf(&outbuf, info_s, strlen(info_s));

                    } else if (q[k].variant_op == p_pct) {
                        asprintf(&info_s, ";GTQ_%u=%f", k,
                                //((float)counts[k][vids[line_idx]])/
                                ((float)counts[k][line_idx])/
                                ((float) id_lens[k]));
                        append_out_buf(&outbuf, info_s, strlen(info_s));
                    } else if (q[k].variant_op == p_maf) {
                        asprintf(&info_s, ";GTQ_%u=%f", k,
                                //((float)counts[k][vids[line_idx]])/
                                ((float)counts[k][line_idx])/
                                (((float) id_lens[k])*2.0));
                        append_out_buf(&outbuf, info_s, strlen(info_s));
                    }

                }
                
	        append_out_buf(&outbuf,"\n",1);
            }
	    bit_i++;
	    if (bit_i == num_fields)
	        break;
        }

        if (bit_i == num_fields)
            break;
    }
    quick_file_delete(&qfile);
    free_out_buf(&outbuf);
}
//}}}

//{{{int query_help()
int query_help()
{
    printf(
"usage:   gqt query -i <wahbm file> \\\n"
"                   [-b <bim file> || -s <bcf file> && -v <vid file>]  \\\n"                    
"                   -c only print number of resulting variants \\\n"
"                   -d <ped database file> \\\n"
"                   -p <population query 1> \\\n"
"                   -g <genotype query 1> \\\n"
"                   -p <population query 2> \\\n"
"                   -g <genotype query 2> \\\n"
"\n"
"A GQT query returns a set of variants that meet some number of population \n"
"and genotype conditions.  Conditions are specified by a population query \n"
"and genotype query pair, where the population query defines the set of\n"
"individuals to consider and the genotype query defines a filter on that\n"
"population.  The result is the set of variants within that sub-population\n"
"that meet the given conditions.  For example, to find the variants that are\n"
"heterozygous in the GBR population the query pair would be:\n"
"\n"
"\t-p \"Population = 'GBR'\" -g \"HET\"\n"
"\n"
"Any number of query pairs can be included, to further refine that set of\n"
"variants.  For example, to find the variants that are heterozygous in at \n"
"least 10 individuals from the GBR population, and are homozygous reference \n"
"in the TSI population the two query pairs would be:\n"
"\n"
"\t-p \"Population = 'GBR'\" -g \"count(HET) >= 10\" \\\n"
"\t-p \"Population = 'GBR'\" -g \"HOMO_REF\"\n"
"\n"
"Population queries are based on the PED file that is associated with the\n"
"genotypes, and any column in that PED file can be part of the query.  For\n"
"example, a PED file that includes the \"Paternal_ID\" and \"Gender\" fields\n"
"(where male = 1 and female = 2) could be queried by:\n"
"\n"
"\t-p \"Paternal_ID = 'NA12878' AND Gender = 2\"\n"
"\n"
"Genotype queries can either be direct genotype filters or count-based \n"
"filters.  To get the variants that are heterozygous in every member of the\n"
"population the query would be:\n"
"\n"
"\t-g \"HET\"\n"
"\n"
"Or to get the variants that are either heterozygous or homozygous alternate\n"
"in every member the query would be:\n"
"\n"
"\t-g \"HET HOMO_ALT\"\n"
"\n"
"Count based filters used the \"count()\" operator that takes a genotype \n"
"list as a parameter followed by some condition.  For example, to find the\n"
"variants that are either heterozygous or homozygous alternate in no more\n"
"than 10 individuals the query would be\n"
"\n"
"\t-g \"count(HET HOMO_ALT) < 10\"\n");
    return 1;
}
//}}}
