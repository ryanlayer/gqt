#define _GNU_SOURCE
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <inttypes.h>

#include "genotq.h"
#include "vid.h"
#include "wahbm_in_place.h"
#include "ped.h"
#include "wahbm.h"
#include "bcf.h"
#include "off.h"
#include "parse_q.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

int query_help(int exit_code);

void print_query_result_offset(uint32_t *mask,
                               uint32_t mask_len,
                               uint32_t *vids,
                               struct gqt_query *q,
                               uint32_t **counts,
                               uint32_t *id_lens,
                               uint32_t *U_R,
                               uint32_t U_R_len,
                               char **id_query_list,
                               char **gt_query_list,
                               uint32_t num_qs,
                               uint32_t num_fields,
                               char *off_file_name,
                               char *source_file,
                               char *full_cmd);

void print_query_result_bim(uint32_t *mask,
                        uint32_t mask_len,
                        uint32_t *vids,
                        struct gqt_query *q,
                        uint32_t **counts,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *bim,
                        char *full_cmd);
int query_cmp(uint32_t value,
              int op_condition,
              int condition_value);


#if 0
void get_bcf_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        struct gqt_query *q,
                        char **id_query_list,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *vid_file_name,
                        char *bcf_file_name,
                        int bcf_output);
#endif

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

//{{{ int query(int argc, char **argv, char *full_cmd)
int query(int argc, char **argv, char *full_cmd)
{
    if (argc < 2) return query_help(EX_USAGE);

    int c;
    char *input_file_name=NULL,
         *gqt_file_name=NULL,
         *ped_db_file_name=NULL,
         *bim_file_name=NULL,
         *off_file_name=NULL,
         *bcf_file_name=NULL,
         *vid_file_name=NULL;
    int c_is_set = 0,
        i_is_set = 0,
        id_q_count = 0,
        gt_q_count = 0,
        d_is_set = 0,
        v_is_set = 0,
        V_is_set = 0,
        G_is_set = 0,
        B_is_set = 0,
        O_is_set = 0,
        S_is_set = 0,
        bcf_output = 0;

    char *id_query_list[100];
    char *gt_query_list[100];

    //{{{ parse cmd line opts
    while ((c = getopt (argc, argv, "chvi:p:g:d:B:V:G:O:S:")) != -1) {
        switch (c) {
        case 'v':
            v_is_set = 1;
            break;
        case 'c':
            c_is_set = 1;
            break;
        case 'i':
            i_is_set = 1;
            input_file_name = optarg;
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
            ped_db_file_name = optarg;
            break;
        case 'B':
            B_is_set = 1;
            bim_file_name = optarg;
            break;
        case 'V':
            V_is_set = 1;
            vid_file_name = optarg;
            break;
        case 'G':
            G_is_set = 1;
            gqt_file_name = optarg;
            break;
        case 'S':
            S_is_set = 1;
            bcf_file_name = optarg;
            break;
        case 'O':
            O_is_set = 1;
            off_file_name = optarg;
            break;
        case 'h':
            return query_help(EX_OK);
        case '?':
            if ( (optopt == 'i') ||
                    (optopt == 'p') ||
                    (optopt == 'g') ||
                    (optopt == 'd') ||
                    (optopt == 'b') )
                fprintf (stderr,
                        "Option -%c requires an argument.\n",
                         optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            return query_help(EX_USAGE);
        default:
            return query_help(EX_OK);
        }
    }

    if (i_is_set == 0) {
        fprintf (stderr, "No input file given (vcf.gz/bcf/gqt).\n");
        return query_help(EX_NOINPUT);
    }

    if (strlen(input_file_name) < 4) {
        fprintf (stderr,
                 "Cannot determine input file type for file '%s'.\n",
                 input_file_name);
        fprintf (stderr,
                 "NOTE: The file must have a supported extension "
                 "(vcf.gz/bcf/gqt).\n");
        return query_help(EX_NOINPUT);
    }
    // See if -i is a gqt file or a bcf file
    char *input_file_type = (char *)malloc(4*sizeof(char));

    strncpy(input_file_type,
            input_file_name + ( strlen(input_file_name) - 3),
            3 * sizeof(char));
    input_file_type[3] = '\0';

    if (strcmp(input_file_type, "gqt") == 0 ) {
        G_is_set = 1;
        gqt_file_name = input_file_name;
    } else if (strcmp(input_file_type, "bcf") == 0 ) {
        S_is_set = 1;
        bcf_file_name = input_file_name;
   } else {
        if (strlen(input_file_name) < 8) {
            fprintf (stderr,
                    "Cannot determine input file type for file '%s'.\n",
                    input_file_name);
            fprintf (stderr,
                     "NOTE: The file must have a supported extension "
                     "(vcf.gz/bcf/gqt).\n");
            return query_help(EX_NOINPUT);
        }

        free(input_file_type);
        input_file_type = (char *)malloc(6*sizeof(char));

        strncpy(input_file_type,
                input_file_name + ( strlen(input_file_name) - 6),
                6 * sizeof(char));

        if (strcmp(input_file_type, "vcf.gz") == 0 ) {
            S_is_set = 1;
            bcf_file_name = input_file_name;
        } else {
            fprintf (stderr,
                    "Cannot determine input file type for file '%s'.\n",
                    input_file_name);
            fprintf (stderr,
                     "NOTE: The file must have a supported extension "
                     "(vcf.gz/bcf/gqt).\n");
            return query_help(EX_NOINPUT);
        }
    }

    // BCF/VCFGZ file is set
    if (S_is_set == 1) {
        if ( access( bcf_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing BCF file '%s'", bcf_file_name);

        // GQT is not, autodetect
        if (G_is_set == 0) {
            gqt_file_name  = (char*)malloc(strlen(bcf_file_name) + 5); 
            if (!gqt_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(gqt_file_name, bcf_file_name);
            strcat(gqt_file_name, ".gqt");

            if ( access( gqt_file_name, F_OK) != -1 ) {
                G_is_set = 1;
            } else {
                fprintf(stderr,
                        "Auto detect failure: GQT file '%s' not found\n",
                        gqt_file_name);
                return query_help(EX_NOINPUT);
            }
        } else {
            if ( access( gqt_file_name, F_OK) == -1 ) {
                fprintf(stderr, "GQT file '%s' not found\n", gqt_file_name);
                return query_help(EX_NOINPUT);
            }
        }

        // PED DB is not, autodetect
        if (d_is_set == 0) {
            ped_db_file_name  = (char*)malloc(strlen(bcf_file_name) + 4); 
            if (!ped_db_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(ped_db_file_name, bcf_file_name);
            strcat(ped_db_file_name, ".db");

            if ( access( ped_db_file_name, F_OK) != -1 ) {
                d_is_set = 1;
            } else {
                fprintf(stderr,
                        "Auto detect failure: DB file '%s' not found\n",
                        ped_db_file_name);
                return query_help(EX_NOINPUT);
            }
        }

        // VID is not, autodetect
        if (V_is_set == 0) {
            vid_file_name  = (char*)malloc(strlen(bcf_file_name) + 5); 
            if (!vid_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(vid_file_name, bcf_file_name);
            strcat(vid_file_name, ".vid");

            if ( access( vid_file_name, F_OK) != -1 ) {
                V_is_set = 1;
            } else {
                fprintf(stderr,
                        "Auto detect failure: VID file '%s' not found\n",
                        vid_file_name);
                return query_help(EX_NOINPUT);
            }
        }

        // Try and find the BIM file, okay if not there (for now)
        if (B_is_set == 0) {
            bim_file_name  = (char*)malloc(strlen(bcf_file_name) + 5); 
            if (!bim_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(bim_file_name, bcf_file_name);
            strcat(bim_file_name, ".bim");
            if ( access( bim_file_name, F_OK) != -1 ) 
                B_is_set = 1;
        }

        // Try and find the OFF file, okay if not there (for now)
        if (O_is_set == 0) {
            off_file_name  = (char*)malloc(strlen(bcf_file_name) + 5); 
            if (!off_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(off_file_name, bcf_file_name);
            strcat(off_file_name, ".off");

            if ( access( off_file_name, F_OK) != -1 ) 
                O_is_set = 1;
        } 
    } else if (G_is_set == 1) { 
        if ( access( gqt_file_name, F_OK) == -1 ) {
            fprintf(stderr, "GQT file '%s' not found\n", gqt_file_name);
            return query_help(EX_NOINPUT);
        }

        // Try and find the BIM file, okay if not there (for now)
        if (B_is_set == 0) {
            bim_file_name  = (char*)
                    malloc((1+strlen(gqt_file_name))*sizeof(char)); 
            if (!bim_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(bim_file_name, gqt_file_name);
            strcpy(bim_file_name + strlen(bim_file_name) - 3, "bim");
            if ( access( bim_file_name, F_OK) != -1 ) 
                B_is_set = 1;
        } 

        if (O_is_set == 0) {
            off_file_name  = (char*)
                malloc((1+strlen(gqt_file_name))*sizeof(char)); 
            if (!off_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(off_file_name, gqt_file_name);
            strcpy(off_file_name + strlen(off_file_name) - 3, "off");
            if ( access( off_file_name, F_OK) != -1 ) 
                O_is_set = 1;
        } 


        if (V_is_set == 0) {
            vid_file_name  = (char*)
                malloc((1+strlen(gqt_file_name))*sizeof(char)); 
            if (!vid_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(vid_file_name, gqt_file_name);
            strcpy(vid_file_name + strlen(vid_file_name) - 3, "vid");

            if ( access( vid_file_name, F_OK) != -1 ) {
                V_is_set = 1;
            } else {
                fprintf(stderr,
                        "Auto detect failure: VID file '%s' not found\n",
                        vid_file_name);
                return query_help(EX_NOINPUT);
            }
        }

        if (d_is_set == 0) {
            ped_db_file_name  = (char*)
                malloc((1+strlen(gqt_file_name))*sizeof(char)); 
            if (!ped_db_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(ped_db_file_name, gqt_file_name);
            strcpy(ped_db_file_name + strlen(gqt_file_name) - 3, "db\0");

            if ( access( ped_db_file_name, F_OK) != -1 ) {
                d_is_set = 1;
            } else {
                fprintf(stderr,
                        "Auto detect failure: PED DB file '%s' not found\n",
                        ped_db_file_name);
                return query_help(EX_NOINPUT);
            }
        }
    }  else {
        fprintf(stderr,
                "Neither GQT or BCF/VCF.GZ file given.\n");
        return query_help(EX_NOINPUT);
    } 

    if (B_is_set == 1) {
        if ( access( bim_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing BIM file '%s'", bim_file_name);
    }

    if (O_is_set == 1) {
        if ( access( off_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing OFF file '%s'", bim_file_name);
    }

    if (V_is_set == 0) {
        fprintf(stderr, "VID file is not set\n");
        return query_help(EX_NOINPUT);
    } else {
        if ( access( vid_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing VID file '%s'", vid_file_name);
    }

    if (d_is_set == 0) {
        fprintf(stderr, "DB file is not set\n");
        return query_help(EX_NOINPUT);
    } else {
        if ( access( ped_db_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing DB file '%s'", ped_db_file_name);
    }

    if (gt_q_count != id_q_count) {
        fprintf(stderr, 
                "Mismatched number of individual and genotype query strings\n");
        return query_help(EX_USAGE);
    }


    if ((B_is_set == 0) && (O_is_set == 0) && (c_is_set == 0)) {
        fprintf(stderr, 
                "Must set either BIM or OFF files when doing anything "
                "other than counting.\n");
        return query_help(EX_USAGE);
    }

    if ( (v_is_set == 1) && ((O_is_set == 0) || (S_is_set == 0)) ) {
        fprintf(stderr, 
                "To get genotypes source BCF/VCF.GZ and OFF files "
                "must be set.\n");
        return query_help(EX_USAGE);
    }

    if ( (v_is_set == 0) && (B_is_set == 0) ) {
        fprintf(stderr, 
                "To get variant data only BIM file must be set.\n");
        return query_help(EX_USAGE);
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

    struct wahbm_file *wf = open_wahbm_file(gqt_file_name);
    struct vid_file *vid_f = open_vid_file(vid_file_name);
    load_vid_data(vid_f);

    //uint32_t num_ints = (wf.num_fields + 32 - 1)/ 32;
    uint32_t num_ints = (wf->gqt_header->num_variants + 32 - 1)/ 32;
    uint32_t len_ints;
    uint32_t *U_R = NULL;
    uint32_t U_R_len = 0;

    for (i = 0; i < gt_q_count; ++i) {
        uint32_t len_count_R;
        uint32_t *R;
        /* 
         * Submit the population query to the PED database and get back both
         * the list of of ids in R and the length of R in id_lens[i]
         */
        id_lens[i] = resolve_ind_query(&R,
                                      id_query_list[i],
                                      ped_db_file_name);

        uint32_t *tmp_U_R = (uint32_t *)
                realloc(U_R, (U_R_len + id_lens[i]) * sizeof(uint32_t));
        if (!tmp_U_R)
            err(EX_OSERR, "malloc error");
        else
            U_R = tmp_U_R;


        for (j = 0; j < id_lens[i]; ++j) {
            U_R[U_R_len] = R[j];
            U_R_len += 1;
        }

        // Enforce that the offsets of the relevant samples is 
        // within the number of samples in the GQT index.
        if (id_lens[i] > wf->gqt_header->num_samples) {
            fprintf(stderr, 
                    "ERROR: there are more samples in the PED database (%d) "
                    "that match this condition \nthan there are in the GQT "
                    "index (%d).  Perhaps your PED file is a superset of "
                    "the\nsamples in your VCF/BCF file?\n", 
                    id_lens[i], 
                    wf->gqt_header->num_samples);
            return 1;
        }

        uint32_t low_v = 0, high_v = 0;

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
                mapped_counts[i][vid_f->vids[j]] = counts[i][j];

            gt_mask[i] = (uint32_t *) malloc(num_ints * sizeof(uint32_t));
            if (!gt_mask[i])
                err(EX_OSERR, "malloc error");

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

    if (U_R == NULL) {
        U_R_len= resolve_ind_query(&U_R,
                                   "",
                                   ped_db_file_name);
    }

    // Get the uniq elements in place
    qsort(U_R, U_R_len, sizeof(uint32_t), compare_uint32_t);
    for (i = j = 0; i < U_R_len; i++)
        if (U_R[i] != U_R[j]) 
            U_R[++j] = U_R[i];
    U_R_len = j + 1;


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

        if (masked_vid_count <= wf->gqt_header->num_variants)
            printf("%u\n", masked_vid_count);
        else
            printf("%u\n", wf->gqt_header->num_variants);

    } else if ( (B_is_set == 1) || (O_is_set == 1)){

        uint32_t *mapped_mask = (uint32_t *) calloc(num_ints,sizeof(uint32_t));

        uint32_t v,p,leading_zeros=32, hit;
        for (i = 0; i < num_ints; ++i) {
            if (final_mask[i] != 0) {
                v = final_mask[i];
                p = popcount(v);
                for (j = 0; j < p; ++j) {
                    leading_zeros = __builtin_clz(v);

                    if (i*32 + leading_zeros + 1 > wf->gqt_header->num_variants)
                        break;

                    hit = vid_f->vids[leading_zeros + i*32];

                    mapped_mask[hit/32] |= 1 << (31-hit%32);
                    v &= ~(1 << (32 - leading_zeros - 1));
                }
            }
            if (i*32 + leading_zeros + 1 > wf->gqt_header->num_variants)
                break;
        }

        if (v_is_set == 0)
            print_query_result_bim(mapped_mask,
                               num_ints,
                               vid_f->vids,
                               q,
                               mapped_counts,
                               id_lens,
                               gt_q_count,
                               wf->gqt_header->num_variants,
                               bim_file_name,
                               full_cmd);
        else
            print_query_result_offset(mapped_mask,
                                      num_ints,
                                      vid_f->vids,
                                      q,
                                      mapped_counts,
                                      id_lens,
                                      U_R,
                                      U_R_len,
                                      id_query_list,
                                      gt_query_list,
                                      gt_q_count,
                                      wf->gqt_header->num_variants,
                                      off_file_name,
                                      bcf_file_name,
                                      full_cmd);
    }

    for (j = 0; j < gt_q_count; ++j) {
        free(gt_mask[j]);
        if ( (q[j].variant_op == p_count) || 
             (q[j].variant_op == p_pct) ||
             (q[j].variant_op == p_maf) )
            free(counts[j]);
    }

    destroy_vid_file(vid_f);
    destroy_wahbm_file(wf);
    return 0;
}
//}}}
 
//{{{ void print_query_result_bim(uint32_t *mask,
void print_query_result_bim(uint32_t *mask,
                            uint32_t mask_len,
                            uint32_t *vids,
                            struct gqt_query *q,
                            uint32_t **counts,
                            uint32_t *id_lens,
                            uint32_t num_qs,
                            uint32_t num_fields,
                            char *bim,
                            char *full_cmd)
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

    char *info_s = NULL;;

    int r = asprintf(&info_s, 
                     "##%s_queryVersion=%s\n"
                     "##%s_queryCommand=%s\n",
                     PROGRAM_NAME, VERSION,
                     PROGRAM_NAME, full_cmd);
    if (r == -1) err(EX_OSERR, "asprintf error");
                      
    append_out_buf(&outbuf, info_s, strlen(info_s));

    for (k=0; k < num_qs; k++) {
        if ( q[k].variant_op == p_count ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Integer,"
                         "Description=\"GQT count result from query "
                         "%u\">\n",
                         k, k);
            if (r == -1) err(EX_OSERR, "asprintf error");
            append_out_buf(&outbuf, info_s, strlen(info_s));
        } else if ( q[k].variant_op == p_pct ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Float,"
                         "Description=\"GQT percent result from query "
                         "%u\">\n",
                         k, k);
            if (r == -1) err(EX_OSERR, "asprintf error");
            append_out_buf(&outbuf, info_s, strlen(info_s));
        } else if ( q[k].variant_op == p_maf ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Float,"
                         "Description=\"GQT maf result from query "
                         "%u\">\n",
                         k, k);
            if (r == -1) err(EX_OSERR, "asprintf error");
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
                        r = asprintf(&info_s,
                                     ";GQT_%u=%u",
                                     k,
                                     counts[k][line_idx]);
                                //counts[k][vids[line_idx]]);
                        if (r == -1) err(EX_OSERR, "asprintf error");
                        append_out_buf(&outbuf, info_s, strlen(info_s));

                    } else if (q[k].variant_op == p_pct) {
                        r = asprintf(&info_s, ";GQT_%u=%f", k,
                                //((float)counts[k][vids[line_idx]])/
                                ((float)counts[k][line_idx])/
                                ((float) id_lens[k]));
                        if (r == -1) err(EX_OSERR, "asprintf error");
                        append_out_buf(&outbuf, info_s, strlen(info_s));
                    } else if (q[k].variant_op == p_maf) {
                        r = asprintf(&info_s, ";GQT_%u=%f", k,
                                //((float)counts[k][vids[line_idx]])/
                                ((float)counts[k][line_idx])/
                                (((float) id_lens[k])*2.0));
                        if (r == -1) err(EX_OSERR, "asprintf error");
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

//{{{ void print_query_result_offset(uint32_t *mask,
void print_query_result_offset(uint32_t *mask,
                               uint32_t mask_len,
                               uint32_t *vids,
                               struct gqt_query *q,
                               uint32_t **counts,
                               uint32_t *id_lens,
                               uint32_t *U_R,
                               uint32_t U_R_len,
                               char **id_query_list,
                               char **gt_query_list,
                               uint32_t num_qs,
                               uint32_t num_fields,
                               char *off_file_name,
                               char *source_file,
                               char *full_cmd)
{
    struct off_file *off_f = open_off_file(off_file_name);
    struct bcf_file bcf_f = init_bcf_file(source_file);

    char *sample_names = NULL;

    uint32_t i,j,k,line_idx,bytes, bit_i = 0;
    int r;
    for (i = 0; i < U_R_len; ++i) {
        if (i == 0 )
            r = asprintf(&sample_names,
                         "%s",
                         bcf_f.hdr->samples[U_R[i]]);
        else
            r = asprintf(&sample_names,
                         "%s,%s",
                         sample_names,
                         bcf_f.hdr->samples[U_R[i]]);
        if (r == -1)
            err(EX_OSERR, "asprintf error");
    }

    if (bcf_hdr_set_samples(bcf_f.hdr, sample_names, 0) != 0)
        errx(EX_DATAERR, "Error setting samples: %s\n", source_file);

    char *info_s;
 
    for (i = 0; i < num_qs; i++) {
        if ( q[i].variant_op == p_count ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Integer,"
                         "Description=\"GQT count result from "
                         "phenotype:'%s' genotype:'%s'\">",
                         i, id_query_list[i], gt_query_list[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");

            if (bcf_hdr_append(bcf_f.hdr, info_s) != 0)
                errx(EX_DATAERR, "Error updating header: %s\n", source_file);

        } else if ( q[i].variant_op == p_pct ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Float,"
                         "Description=\"GQT percent result from "
                         "phenotype:'%s' genotype:'%s'\">",
                         i, id_query_list[i], gt_query_list[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");

            if (bcf_hdr_append(bcf_f.hdr, info_s) != 0)
                errx(EX_DATAERR, "Error updating header: %s\n", source_file);

        } else if ( q[i].variant_op == p_maf ) {
            r = asprintf(&info_s, "##INFO=<ID=GQT_%u,Number=1,Type=Float,"
                         "Description=\"GQT maf result from "
                         "phenotype:'%s' genotype:'%s'\">",
                         i, id_query_list[i], gt_query_list[i]);

            if (bcf_hdr_append(bcf_f.hdr, info_s) != 0)
                errx(EX_DATAERR, "Error updating header: %s\n", source_file);
        }

    }

    r = asprintf(&info_s, "##%s_queryVersion=%s", PROGRAM_NAME, VERSION);
    if (r == -1) err(EX_OSERR, "asprintf error");

    if (bcf_hdr_append(bcf_f.hdr, info_s) != 0)
        errx(EX_DATAERR, "Error updating header: %s\n", source_file);

    r = asprintf(&info_s, "##%s_queryCommand=%s", PROGRAM_NAME, full_cmd);
    if (r == -1) err(EX_OSERR, "asprintf error");

    if (bcf_hdr_append(bcf_f.hdr, info_s) != 0)
        errx(EX_DATAERR, "Error updating header: %s\n", source_file);

    htsFile *out_f = hts_open("-","w");
    if ( !out_f )
        err(EX_DATAERR, "Could open output file");

    bcf_hdr_write(out_f, bcf_f.hdr);

    bcf_f.line = bcf_init1();

    for (i=0; i < mask_len; ++i) {
        bytes = mask[i];
	if (bytes == 0)
            continue; /* skip a bunch of ops if you can */
        for (j=0; j < 32; j++) {
            if (bytes & 1 << (31 - j)) {
	        line_idx = i*32+j;

                r = goto_bcf_line(&bcf_f, off_f, line_idx);

                if (r == -1) 
                    err(EX_NOINPUT,
                        "Error seeking file '%s'", bcf_f.file_name);

                r = get_bcf_line(&bcf_f);
                if (r == -1) 
                    err(EX_NOINPUT,
                        "Error reading file '%s'", bcf_f.file_name);

                for (k=0; k < num_qs; k++) {
                    r = asprintf(&info_s, "GQT_%u", k);
                    if (r == -1)
                        err(EX_OSERR, "asprintf error");

                    if ( q[k].variant_op == p_count ) {
                        int32_t v = counts[k][line_idx];
                        if (bcf_update_info_int32(bcf_f.hdr,
                                                  bcf_f.line,
                                                  info_s,
                                                  &v,
                                                  1) != 0)
                            errx(EX_DATAERR,
                                 "Error adding to info field: %s\n",
                                 bcf_f.file_name);
                    } else if (q[k].variant_op == p_pct) {
                        float v = ((float)counts[k][line_idx])/
                                    ((float) id_lens[k]);
                        if (bcf_update_info_float(bcf_f.hdr,
                                                  bcf_f.line,
                                                  info_s,
                                                  &v,
                                                  1) != 0)
                            errx(EX_DATAERR,
                                 "Error adding to info field: %s\n",
                                 bcf_f.file_name);

                    } else if (q[k].variant_op == p_maf) {
                        float v = ((float)counts[k][line_idx])/
                                    (((float) id_lens[k])*2.0);
                        if (bcf_update_info_float(bcf_f.hdr,
                                                  bcf_f.line,
                                                  info_s,
                                                  &v,
                                                  1) != 0)
                            errx(EX_DATAERR,
                                 "Error adding to info field: %s\n",
                                 bcf_f.file_name);
                    }
                }

                bcf_write(out_f, bcf_f.hdr, bcf_f.line);

            }
	    bit_i++;
	    if (bit_i == num_fields)
	        break;
        }

        if (bit_i == num_fields)
            break;
    }
    hts_close(out_f);
    destroy_off_file(off_f);
}
//}}}

//{{{int query_help()
int query_help(int exit_code)
{
    fprintf(stderr, 
"%s v%s\n"
"usage:   gqt query -i <bcf/vcf or gqt file> \\\n"
"                   -d <ped database file> \\\n"
"                   -c only print number of resulting variants \\\n"
"                   -v print genotypes (from the source bcf/vcf)\\\n"
"                   -B <bim file> (opt.)\\\n"
"                   -O <off file> (opt.)\\\n"
"                   -V <vid file> (opt.)\\\n"
"                   -G <gqt file> (opt.)\\\n"
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
"\t-p \"Population = 'GBR'\" -g \"HOM_REF\"\n"
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
"\t-g \"HET HOM_ALT\"\n"
"\n"
"Count based filters used the \"count()\" operator that takes a genotype \n"
"list as a parameter followed by some condition.  For example, to find the\n"
"variants that are either heterozygous or homozygous alternate in no more\n"
"than 10 individuals the query would be\n"
"\n"
"\t-g \"count(HET HOM_ALT) < 10\"\n",
PROGRAM_NAME, VERSION);
    return exit_code;
}
//}}}

#if 0 
//{{{ void get_bcf_query_result(uint32_t *mask,
void get_bcf_query_result(uint32_t *mask,
                        uint32_t mask_len,
                        struct gqt_query *q,
                        char **id_query_list,
                        uint32_t *id_lens,
                        uint32_t num_qs,
                        uint32_t num_fields,
                        char *vid_file_name,
                        char *bcf_file_name,
                        int bcf_output)
{

    /* The VID file contains the line numbers of the variants after they have
     * been sorted.  To reach back into the BCF file to print the metadata
     * associated with the variants marked in the mask, we need to create a
     * sorted list of line numbers we want.  So first we intersect the VID file
     * and the mask, then sort it.
     */
    /*
    FILE *vid_f = fopen(vid_file_name, "rb");
    if (!vid_f)
        err(EX_NOINPUT, "Cannot read file\"%s\"", vid_file_name);

    uint32_t *vids = (uint32_t *) malloc(num_fields*sizeof(uint32_t));
    if (!vids )
        err(EX_OSERR, "malloc error");

    size_t fr = fread(vids, sizeof(uint32_t), num_fields, vid_f);
    check_file_read(vid_file_name, vid_f, num_fields, fr);

    fclose(vid_f);
    */
    struct vid_file *vid_f = open_vid_file(vid_file_name);
    load_vid_data(vid_f);

    uint32_t i, j, masked_vid_count = 0;

    for (i = 0; i < mask_len; ++i)
        masked_vid_count += popcount(mask[i]);

    uint32_t *masked_vids = (uint32_t *)
            malloc(masked_vid_count*sizeof(uint32_t));
    if (!masked_vids )
        err(EX_OSERR, "malloc error");
    uint32_t masked_vid_i = 0;

    for (i = 0; i < mask_len; ++i) {
        uint32_t bytes = mask[i];
	if (bytes == 0)
            continue; /* skip a bunch of ops if you can */
        for (j = 0; j < 32; j++) {
            if (bytes & (1 << (31 - j))) {
                masked_vids[masked_vid_i] = vid_f->vids[i*32 + j];
                masked_vid_i+=1;
            }
        }
        if (masked_vid_i == masked_vid_count)
            break;
    }

    destroy_vid_file(vid_f);

    qsort(masked_vids, masked_vid_count, sizeof(uint32_t), compare_uint32_t);

    htsFile *fp    = hts_open(bcf_file_name,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *line    = bcf_init1();
    //bcf_hdr_set_samples(hdr, print_name_csv, 0);

    htsFile *out;
    if (!bcf_output)
        out = hts_open("-", "w");
    else
        out = hts_open("-", "wb");

    int r = bcf_hdr_write(out, hdr);

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
#endif
