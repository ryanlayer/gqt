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

int ihs_help();

//{{{ int ihs(int argc, char **argv)
int ihs(int argc, char **argv)
{
#if 0
    if (argc < 2) return ihs_help();

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


    //{{{ parse cmd line opts
    while ((c = getopt (argc, argv, "chi:p:g:d:b:v:s:B")) != -1) {
        switch (c) {
        case 'i':
            i_is_set = 1;
            wahbm_file_name = optarg;
            break;
        case 'p':
            id_query_list[id_q_count] = optarg;
            id_q_count += 1;
            break;
        case 'd':
            d_is_set = 1;
            db_file_name = optarg;
            break;
        case 'h':
            return ihs_help();
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
            return ihs_help();
        }
    }

    if (i_is_set == 0) {
        printf("WAHBM file is not set\n");
        return ihs_help();
    }

    if (d_is_set == 0) {
        printf("PED database file is not set\n");
        return ihs_help();
    }

    //}}}

    struct gqt_query q[100];
    uint32_t *gt_mask[100];
    uint32_t *counts[100];
    uint32_t id_lens[100];
    char *id_query_list[100];

    int r, i, j, k;

    struct wah_file wf = init_wahbm_file(wahbm_file_name);
    uint32_t num_ints = (wf.num_fields + 32 - 1)/ 32;
    uint32_t len_ints;

    for (i = 0; i < id_q_count; ++i) {
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

        /*
         * q holds the parameters of each query, first determin the range of 
         * bitmaps to pull
         */
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

        /* User asks for a count of percent */
        if ( ( q[i].variant_op == p_count ) || ( q[i].variant_op == p_pct ) ) {


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


            gt_mask[i] = (uint32_t *) malloc(num_ints * sizeof(uint32_t));

            /* User specifies a condition */
            if ( q[i].op_condition != -1) { 

                /* Since we only find counts, when the user asks for a
                 * perecent, just convert that back to the count that meets the
                 * percent condition
                 */
                float condition_value = q[i].condition_value;
                if ( q[i].variant_op == p_pct )
                    condition_value *= id_lens[i];


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

        /*
        fprintf(stderr,"mask[%u]:", j);
        for (j = 0; j < num_ints; ++j)
            fprintf(stderr,"\t%u", gt_mask[i][j]);
        fprintf(stderr,"\n");
        */
                

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
        print_query_result(final_mask,
                        num_ints,
                        q,
                        counts,
                        id_lens,
                        gt_q_count,
                        wf.num_fields,
                        bim_file_name);
    }

    for (j = 0; j < gt_q_count; ++j) {
        free(gt_mask[j]);
        if ( ( q[j].variant_op == p_count ) || ( q[j].variant_op == p_pct ) )
            free(counts[j]);
    }

    fclose(wf.file);
#endif
    return 0;
}
//}}}

//{{{int query_help()
int ihs_help()
{
    printf(
"usage:   gqt ihs -i <gqt file> \\\n"
"                   -c only print number of resulting variants \\\n"
"                   -d <ped database file> \\\n"
"                   -p <population query 1> \\\n"
"                   -p <population query 2> \\\n"
"\n");
    return 1;
}
//}}}
