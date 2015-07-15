#include <immintrin.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <inttypes.h>
#include <time.h>
#include "genotq.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

void permute(uint32_t *R, uint32_t len_R);

int pop_help(char *op);

uint32_t pca_shared(struct wah_file *wf,
                    char **id_query_list,
                    char *label_field_name,
                    char *db_file_name,
                    char *id_out_file);

uint32_t fst(struct wah_file *wf,
             char **id_query_list,
             uint32_t id_q_count,
             uint32_t *vids,
             char *db_file_name,
             double **mapped_fst);

uint32_t gst(struct wah_file *wf,
             char **id_query_list,
             uint32_t id_q_count,
             uint32_t *vids,
             char *db_file_name,
             double **mapped_gst);

uint32_t calpha(struct wah_file *wf,
                char **id_query_list,
                uint32_t id_q_count,
                uint32_t *vids,
                char *db_file_name,
                uint32_t **case_ctrl_counts,
                uint32_t *n_cases,
                uint32_t *n_ctrls,
                uint32_t N);

void print_calpha_result(uint32_t *R,
                         uint32_t n_cases,
                         uint32_t n_ctrls,
                         uint32_t N,
                         uint32_t num_variants,
                         char *bim);

void print_pop_result(char *op,
                      double *R,
                      uint32_t num_variants,
                      char *bim);

//{{{ int pop(char *op, int argc, char **argv)
int pop(char *op, int argc, char **argv)
{
    if (argc < 2) return pop_help(op);

    int c;
    char *wahbm_file_name=NULL,
         *id_query=NULL,
         *db_file_name=NULL,
         *bim_file_name=NULL,
         *label_file_name=NULL,
         *label_field_name=NULL,
         *vid_file_name=NULL;
    int i_is_set = 0,
        id_q_count = 0,
        d_is_set = 0,
        v_is_set = 0,
        f_is_set = 0,
        l_is_set = 0,
        b_is_set = 0;
    uint32_t N;

    char *id_query_list[100];

    //{{{ parse cmd line opts
    while ((c = getopt (argc, argv, "hi:p:d:b:v:f:l:")) != -1) {
        switch (c) {
        case 'l':
            l_is_set = 1;
            label_file_name = optarg;
            break;
        case 'f':
            f_is_set = 1;
            label_field_name = optarg;
            break;
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
        case 'b':
            b_is_set = 1;
            bim_file_name = optarg;
            break;
        case 'v':
            v_is_set = 1;
            vid_file_name = optarg;
            break;
        case 'h':
            return pop_help(op);
        case '?':
            if ( (optopt == 'i') ||
                 (optopt == 'l') ||
                 (optopt == 'p') ||
                 (optopt == 'f') ||
                 (optopt == 'd') ||
                 (optopt == 'b') )
                fprintf (stderr, "Option -%c requires an argument.\n",
                         optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        default:
            return pop_help(op);
        }
    }

    if (i_is_set == 0) {
        fprintf(stderr, "GQT file is not set\n");
        return pop_help(op);
    } else {
        if ( access( wahbm_file_name, F_OK) == -1 )
            err(EX_NOINPUT, "Error accessing GQT file \"%s\"", wahbm_file_name);
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
            fprintf(stderr,
                   "Auto detect failure: BIM file %s not found\n",
                   bim_file_name);
            return pop_help(op);
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
            fprintf(stderr,
                    "Auto detect failure: VID file %s not found\n",
                    vid_file_name);
            return pop_help(op);
        }
    }

    if ( (i_is_set == 1) && (d_is_set == 0)) {

        int auto_db_file_name_size = asprintf(&db_file_name,
                                              "%s",
                                              wahbm_file_name);
        strcpy(db_file_name + strlen(db_file_name) - 3, "db\0");

        if ( access( db_file_name, F_OK) != -1 ) {
            d_is_set = 1;
        } else {
            fprintf(stderr, "Auto detect failure: PED DB file %s not found\n",
                   db_file_name);
            return pop_help(op);
        }
    }

    if (i_is_set == 0) {
        fprintf(stderr, "GQT file is not set\n");
        return pop_help(op);
    }

    if (v_is_set == 0) {
        fprintf(stderr, "VID file is not set\n");
        return pop_help(op);
    }

    if (b_is_set == 0) {
        fprintf(stderr, "BIM file is not set\n");
        return pop_help(op);
    }

    if (d_is_set == 0) {
        fprintf(stderr, "PED database file is not set\n");
        return pop_help(op);
    }
    //}}}

    // open WAH/GQT file
    struct wah_file wf = init_wahbm_file(wahbm_file_name);

    // open VID file
    FILE *vid_f = fopen(vid_file_name, "rb");
    if (!vid_f)
        err(EX_NOINPUT, "Cannot read file\"%s\"", vid_file_name);

    uint32_t *vids = (uint32_t *) malloc(wf.num_fields*sizeof(uint32_t));
    int r = fread(vids, sizeof(uint32_t), wf.num_fields, vid_f);
    fclose(vid_f);

    uint32_t num_ints = (wf.num_fields + 32 - 1)/ 32;
    uint32_t len_ints;

    uint32_t len_R;

    if (strcmp("pca-shared",op) == 0) {
        if (f_is_set == 0) {
            fprintf(stderr, "Label database field name not set.\n");
            return pop_help(op);
        }

        if (l_is_set == 0) {
            fprintf(stderr, "Label output file name not set.\n");
            return pop_help(op);
        }

        int r = pca_shared(&wf,
                           id_query_list,
                           label_field_name,
                           db_file_name,
                           label_file_name);
    }else if (strcmp("fst",op) == 0) {
        double *R;
        len_R = fst(&wf, id_query_list, id_q_count, vids, db_file_name, &R);
        print_pop_result(op, R, wf.num_fields, bim_file_name);
    } else if (strcmp("gst",op) == 0) {
        double *R;
        len_R = gst(&wf, id_query_list, id_q_count, vids, db_file_name, &R);
        print_pop_result(op, R, wf.num_fields, bim_file_name);
    } else if (strcmp("calpha",op) == 0) {
        /*
        if (N_is_set == 0) {
            printf("Number of permutations is not set\n");
            return pop_help(op);
        }
        */

        uint32_t n_cases, n_ctrls;
        // The first two arrays hold the observed case/control counts, and
        // the remaining hold the permutation results
        uint32_t *R;

        len_R = calpha(&wf,
                       id_query_list,
                       id_q_count, vids,
                       db_file_name, 
                       &R,
                       &n_cases,
                       &n_ctrls,
                       0);

        print_calpha_result(R,
                            n_cases,
                            n_ctrls,
                            0,
                            wf.num_fields,
                            bim_file_name);
        free(R);

    } else {
        fprintf(stderr, "ERROR: Unknown opperation %s", op);
        pop_help("fst|gst");
    }


    fclose(wf.file);
    return 0;
}
//}}}

//{{{int gst(uint32_t id_q_count,
uint32_t gst(struct wah_file *wf,
             char **id_query_list,
             uint32_t id_q_count,
             uint32_t *vids,
             char *db_file_name,
             double **mapped_gst) {

    uint32_t i,j;
    uint32_t id_lens[id_q_count];
    uint32_t *sums[id_q_count];

    uint32_t num_variants = wf->num_fields;
    uint32_t num_samples = wf->num_records;

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
        if (id_lens[i] > wf->num_records) {
            fprintf(stderr, 
                    "ERROR: there are more samples in the PED database (%d) "
                    "that match this condition \nthan there are in the GQT "
                    "index (%d).  Perhaps your PED file is a superset of "
                    "the\nsamples in your VCF/BCF file?\n", 
                    id_lens[i], 
                    wf->num_records);
            return 1;
        }

        uint32_t low_v, high_v;
        low_v = 1;
        high_v = 3;

#ifdef __AVX2__
        len_count_R = avx_sum_range_records_in_place_wahbm(*wf,
                                                           R,
                                                           id_lens[i],
                                                           low_v,
                                                           high_v,
                                                           &(sums[i]));
#else
        len_count_R = sum_range_records_in_place_wahbm(*wf,
                                                       R,
                                                       id_lens[i],
                                                       low_v,
                                                       high_v,
                                                       &(sums[i]));
#endif
        free(R);
    }

    double *Hs = (double *) calloc(num_variants, sizeof(double));
    double *Ht_a = (double *) calloc(num_variants, sizeof(double));
    double *Ht_b = (double *) calloc(num_variants, sizeof(double));
    double *Fst = (double *) calloc(num_variants, sizeof(double));

    for (i = 0; i < id_q_count; ++i) {
        double maf, pop_size = 2*(double)id_lens[i];
        for (j = 0; j < num_variants; ++j) {
            maf = ((double)sums[i][j])/pop_size;
            Hs[j] += 2*maf*(1-maf);
            Ht_a[j] += maf;
            Ht_b[j] += 1-maf;
        }            
    }

    double num_subpops = (double)id_q_count;
    for (i = 0; i < num_variants; ++i) {
        double Ht_i = 2* (Ht_a[i]/num_subpops)*(Ht_b[i]/num_subpops);
        //Fst[i] = (Hs[i])/num_subpops;
        //Fst[i] = Ht_i;
        Fst[i] = (Ht_i - (Hs[i]/num_subpops)) / Ht_i;
    }

    /* Since the variants are in allele freq order, we need to copy
     * the resulting value to an array that is back in the original
     * order
     */
    *mapped_gst = (double *)calloc(num_variants, sizeof(double));
    for ( i = 0; i < num_variants; ++i)
        (*mapped_gst)[vids[i]] = Fst[i];

    free(Hs);
    free(Ht_a);
    free(Ht_b);
    free(Fst);
    for (i = 0; i < id_q_count; ++i)
        free(sums[i]);

    return num_variants;
}
//}}}

//{{{ uint32_t calpha(struct wah_file *wf,
/*
 * case_ctrl_counts is an array of arrays the following values:
 * there is one array per variants
 * each of those arrays has the vaules:
 * 2: observed # of variants in cases
 * 3: observed # of variants in controls
 * 4: observed # of variants in cases
 * 5: observed # of variants in controls
 * ...
 *
 */
uint32_t calpha(struct wah_file *wf,
                char **id_query_list,
                uint32_t id_q_count,
                uint32_t *vids,
                char *db_file_name,
                uint32_t **mapped_case_ctrl_counts,
                uint32_t *n_cases,
                uint32_t *n_ctrls,
                uint32_t N) {

    uint32_t **case_ctrl_counts = 
            (uint32_t **)malloc((N*2 + 2)* sizeof(uint32_t *));

    uint32_t i,j;
    uint32_t *sums[id_q_count];

    uint32_t num_variants = wf->num_fields;
    uint32_t num_samples = wf->num_records;


    uint32_t *case_ids;
    uint32_t num_cases = resolve_ind_query(&case_ids,
                                           id_query_list[0],
                                           db_file_name);

    uint32_t *ctrl_ids;
    uint32_t num_ctrls = resolve_ind_query(&ctrl_ids,
                                           id_query_list[1],
                                           db_file_name);

    uint32_t num_all = num_cases + num_ctrls;
    uint32_t *all_ids = (uint32_t *)malloc(num_all * sizeof(uint32_t));

    for (i = 0; i < num_cases; ++i) 
        all_ids[i] = case_ids[i];

    for (i = 0; i < num_ctrls; ++i)
        all_ids[i+num_cases] = ctrl_ids[i];

    uint32_t low_v, high_v;
    low_v = 1;
    high_v = 3;

    uint32_t *sum_cases;

#ifdef __AVX2__
    uint32_t len_sum_cases = avx_sum_range_records_in_place_wahbm(*wf,
                                                                  case_ids,
                                                                  num_cases,
                                                                  low_v,
                                                                  high_v,
                                                                  &sum_cases);
#else
    uint32_t len_sum_cases = sum_range_records_in_place_wahbm(*wf,
                                                             case_ids,
                                                             num_cases,
                                                             low_v,
                                                             high_v,
                                                             &sum_cases);
#endif
    case_ctrl_counts[0] = sum_cases;

    uint32_t *sum_ctrls;

#ifdef __AVX2__
    uint32_t len_sum_ctrls = avx_sum_range_records_in_place_wahbm(*wf,
                                                                  ctrl_ids,
                                                                  num_ctrls,
                                                                  low_v,
                                                                  high_v,
                                                                  &sum_ctrls);
#else
    uint32_t len_sum_ctrls = sum_range_records_in_place_wahbm(*wf,
                                                              ctrl_ids,
                                                              num_ctrls,
                                                              low_v,
                                                              high_v,
                                                              &sum_ctrls);
#endif


    case_ctrl_counts[1] = sum_ctrls;

    //srand(time(NULL));
    srand(1);
    uint32_t *P_all_ids = (uint32_t *)malloc(num_all*sizeof(uint32_t));
    memcpy(P_all_ids, all_ids, num_all*sizeof(uint32_t));
    uint32_t *P_case_ids, *P_ctrl_ids;    
    P_case_ids = P_all_ids;
    P_ctrl_ids = P_all_ids + num_cases;

    uint32_t ccc_i = 2;

    for ( i = 0; i < N; ++i) {
        permute(P_all_ids, num_all);

        uint32_t *P_sum_cases;

#ifdef __AVX2__
        uint32_t P_len_sum_cases = 
                avx_sum_range_records_in_place_wahbm(*wf,
                                                P_case_ids,
                                                num_cases,
                                                low_v,
                                                high_v,
                                                &P_sum_cases);
#else
        uint32_t P_len_sum_cases = 
                sum_range_records_in_place_wahbm(*wf,
                                                P_case_ids,
                                                num_cases,
                                                low_v,
                                                high_v,
                                                &P_sum_cases);
#endif

        case_ctrl_counts[ccc_i] = P_sum_cases;
        ccc_i += 1;

        uint32_t *P_sum_ctrls;

#ifdef __AVX2__
        uint32_t P_len_sum_ctrls = 
                avx_sum_range_records_in_place_wahbm(*wf,
                                                     P_ctrl_ids,
                                                     num_ctrls,
                                                     low_v,
                                                     high_v,
                                                     &P_sum_ctrls);
#else
        uint32_t P_len_sum_ctrls = 
                sum_range_records_in_place_wahbm(*wf,
                                                 P_ctrl_ids,
                                                 num_ctrls,
                                                 low_v,
                                                 high_v,
                                                 &P_sum_ctrls);
#endif

        case_ctrl_counts[ccc_i] = P_sum_ctrls;
        ccc_i += 1;
    }


    /*
    uint32_t **mapped_case_ctrl_counts = 
            (uint32_t **)malloc((N*2 + 2)* sizeof(uint32_t *));
    */

    /* 
     * at this point we are going to realign the data so that 
     * the values associate with each variant are grouped together
     * mapped_case_ctrl_counts will have num_variants rows
     * and each row is (N*2 + 2) widde
     *
     * case_ctrl_counts has (N*2 + 2) arrarys, each num_variants long
     * 
     * case_ctrl_counts[i][j]
     * will refere to the ith observation for variant j 
     * which corresponds to 
     * mapped_case_ctrl_counts[j*(N*2 + 2) + i]
     *
     */
    *mapped_case_ctrl_counts = 
            (uint32_t *)malloc((num_variants)*(N*2 + 2)* sizeof(uint32_t *));

    for ( i = 0; i < (N*2 + 2); ++i) {
        for ( j = 0; j < num_variants; ++j) {

            // i * (N*2+2) will get us to the right row
            // and vids[j] to the right col
            (*mapped_case_ctrl_counts)[vids[j]*(N*2+2) + i] =
                case_ctrl_counts[i][j];
        }

        free(case_ctrl_counts[i]);
    }
    free(case_ctrl_counts);

    *n_cases = num_cases;
    *n_ctrls = num_ctrls;

    /*
    for (i = 0; i < num_variants; ++i) {
        for (j = 0; j < (N*2 + 2); ++j)
            fprintf(stderr,
                    "%u\t",
                    (*mapped_case_ctrl_counts)[i*(N*2 +2) + j]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    */

    return num_variants;
}
//}}}

//{{{void permute(uint32_t *R, uint32_t len_R)
void permute(uint32_t *R, uint32_t len_R)
{
    uint32_t k;
    for (k = len_R-1; k > 0; k--) {
        int j = rand() % (k+1);
        int temp = R[j];
        R[j] = R[k];
        R[k] = temp;
    }
}
//}}}

//{{{ int fst(struct wah_file *wf,
uint32_t fst(struct wah_file *wf,
             char **id_query_list,
             uint32_t id_q_count,
             uint32_t *vids,
             char *db_file_name,
             double **mapped_fst) {

    uint32_t i,j;
    uint32_t id_lens[id_q_count];
    uint32_t *sums[id_q_count];
    uint32_t *counts[id_q_count];

    uint32_t num_variants = wf->num_fields;
    uint32_t num_samples = wf->num_records;

    for (i = 0; i < id_q_count; ++i) {
        uint32_t len_count_R, len_sum_R;
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
        if (id_lens[i] > wf->num_records) {
            fprintf(stderr, 
                    "ERROR: there are more samples in the PED database (%d) "
                    "that match this condition \nthan there are in the GQT "
                    "index (%d).  Perhaps your PED file is a superset of "
                    "the\nsamples in your VCF/BCF file?\n", 
                    id_lens[i], 
                    wf->num_records);
            return 1;
        }

        uint32_t low_v, high_v;
        low_v = 1;
        high_v = 3;

#ifdef __AVX2__
        len_sum_R = avx_sum_range_records_in_place_wahbm(*wf,
                                                         R,
                                                         id_lens[i],
                                                         low_v,
                                                         high_v,
                                                         &(sums[i]));
#else
        len_sum_R = sum_range_records_in_place_wahbm(*wf,
                                                     R,
                                                     id_lens[i],
                                                     low_v,
                                                     high_v,
                                                     &(sums[i]));
#endif
        low_v = 1;
        high_v = 2;

#ifdef __AVX2__
        len_count_R = avx_count_range_records_in_place_wahbm(*wf,
                                                             R,
                                                             id_lens[i],
                                                             low_v,
                                                             high_v,
                                                             &(counts[i]));
#else
        len_count_R = count_range_records_in_place_wahbm(*wf,
                                                         R,
                                                         id_lens[i],
                                                         low_v,
                                                         high_v,
                                                         &(counts[i]));
#endif
    }

    //double r = 2.0*(double) num_samples;
    double r = (double) id_q_count;

    // Average sample size
    double nbar = 0;
    for (i = 0; i < id_q_count; ++i) 
        //nbar += (2.0 * (double)id_lens[i])/r;
        nbar += ((double)id_lens[i])/ r;

    double nc = 0;

    for (i = 0; i < id_q_count; ++i) 
        nc += (((double)id_lens[i])*((double)id_lens[i]))/(r*nbar);
    nc = r*nbar - nc;
    nc /= (r - 1.0);

    /* 
     * pbar is the average sample freq of
     *
     * pbar = SUM ( (n_i * p_ij) / (r * nbar) )
     *
     * counts gives us the number over alt alleles for each subpop
     * if id_lens[i] is the number of indv in pop i then
     * let p_ij bet the allele freq of pop i at loc j which is given by
     *      p_i = counts[i][j] / (2.0*id_lens[i])
     */
    double *pbar = (double *) calloc(num_variants, sizeof(double));
    double p_ij, n_i, h_ij;
    for (j = 0; j < num_variants; ++j) {
        for (i = 0; i < id_q_count; ++i) {
            n_i = ((double)id_lens[i]);
            p_ij = ((double)sums[i][j])/(2.0* n_i);
            pbar[j] += (((double)id_lens[i])*p_ij)/(r * nbar);
        }
    }

    /*
     * the sample variant of allele over populations
     * and
     * the average heterozygote freq for the allele
     */
    double *ssqr = (double *) calloc(num_variants, sizeof(double));
    double *hbar = (double *) calloc(num_variants, sizeof(double));
    double *a = (double *) calloc(num_variants, sizeof(double));
    double *b = (double *) calloc(num_variants, sizeof(double));
    double *c = (double *) calloc(num_variants, sizeof(double));
    double *fst = (double *) calloc(num_variants, sizeof(double));

    for (j = 0; j < num_variants; ++j) {
        for (i = 0; i < id_q_count; ++i) {
            n_i = ((double)id_lens[i]);
            p_ij = ((double)sums[i][j])/(2.0*n_i);
            h_ij = ((double)counts[i][j])/n_i;

            ssqr[j] += (n_i * (p_ij - pbar[j]) * (p_ij - pbar[j])) /
                       ((r - 1.0)*nbar);

            hbar[j] += (n_i * h_ij) / ( r * nbar);
        }
        // My understanding of the W-C equations
        a[j] =  (nbar/nc) *
                (   ssqr[j] - (1/(nbar-1.0)) *
                    (   pbar[j]*(1.0 - pbar[j]) -
                        ((r - 1.0)/r)*ssqr[j] - 
                        (0.25)*hbar[j]
                    )
                );
        b[j] =  (nbar/(nbar-1.0)) *
                (   pbar[j]*(1.0 - pbar[j]) -
                    ((r - 1.0)/r)*ssqr[j] -
                    ((2.0*nbar - 1.0)/(4.0*nbar))*hbar[j]
                );
        c[j] = 0.5 * hbar[j];

        fst[j] = a[j]/(a[j] + b[j] + c[j]);
    }

    *mapped_fst = (double *)calloc(num_variants, sizeof(double));
    for ( i = 0; i < num_variants; ++i)
        (*mapped_fst)[vids[i]] = fst[i];

    free(pbar);
    free(ssqr);
    free(hbar);
    free(a);
    free(b);
    free(c);
    free(fst);

    for (i = 0; i < id_q_count; ++i) {
        free(sums[i]);
        free(counts[i]);
    }
    return num_variants;
}
//}}}


//{{{ uint32_t shared(struct wah_file *wf,
uint32_t pca_shared(struct wah_file *wf,
                    char **id_query_list,
                    char *label_field_name,
                    char *db_file_name,
                    char *id_out_file)
{
    uint32_t *record_ids;
    uint32_t num_records = resolve_ind_query(&record_ids,
                                             id_query_list[0],
                                             db_file_name);

    char **labels;
    uint32_t num_labels = resolve_label_query(&labels,
                                              label_field_name,
                                             id_query_list[0],
                                             db_file_name);

    FILE *f = fopen(id_out_file, "w");
    if (!f)
        err(EX_CANTCREAT, "Cannot write to file\"%s\"", id_out_file);

    uint32_t i;
    for (i = 0; i < num_labels; ++i){
        fprintf(f, "%s\n", labels[i]);
        free(labels[i]);
    }

    fclose(f);
 
    uint32_t r = wahbm_shared_by_name_subpop(wf,
                                             record_ids,
                                             num_records);
    return 0;
}
//}}}

//{{{ void print_pop_result(uint32_t *mask,
void print_pop_result(char *op,
                      double *R,
                      uint32_t num_variants,
                      char *bim)
{
    uint32_t i;

    struct quick_file_info qfile;
    struct output_buffer outbuf;
    char pct[50];

    init_out_buf(&outbuf, NULL);
    quick_file_init(bim, &qfile);

    append_out_buf(&outbuf,
                   qfile.main_buf,
                   qfile.header_len);

    char *info_s;

    asprintf(&info_s,
             "##INFO=<ID=GQT_%s,Number=1,Type=Float,"
             "Description=\"GQT %s\">\n",
             op,op);
    append_out_buf(&outbuf, info_s, strlen(info_s));

    char last_header_line[]="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    append_out_buf(&outbuf, last_header_line, strlen(last_header_line));

    for (i=0; i < num_variants; ++i) {
        append_out_buf(&outbuf,
                       qfile.lines[i],
                       qfile.line_lens[i]-1);
        asprintf(&info_s,
                 ";GQT_%s=%f",
                 op,
                 R[i]);
        append_out_buf(&outbuf, info_s, strlen(info_s));
	append_out_buf(&outbuf,"\n",1);
    }

    quick_file_delete(&qfile);
    free_out_buf(&outbuf);
}
//}}}

//{{{ void print_calpha_result(uint32_t *R,
void print_calpha_result(uint32_t *R,
                         uint32_t n_cases,
                         uint32_t n_ctrls,
                         uint32_t N,
                         uint32_t num_variants,
                         char *bim)
{



    uint32_t i,j;
    uint32_t max_int = 0;

    // To reuse the same char* for each of the variants, we need to make sure
    // to allocate enough space to hold all the digits and commas 
    /*
    for (i = 0; i < (num_variants * (N*2 +2)); ++i) {
        if (max_int < R[i])
            max_int = R[i];
    }

    uint32_t max_int_len = (max_int ==0) ? 1 : (int)log10(max_int) + 1;

    uint32_t P_csv_len = (max_int_len + N*2);
    char *P_csv = (char *)malloc(P_csv_len*sizeof(char));
    uint32_t P_csv_i;
    */

    int l;

    struct quick_file_info qfile;
    struct output_buffer outbuf;
    char pct[50];

    init_out_buf(&outbuf, NULL);
    quick_file_init(bim, &qfile);

    append_out_buf(&outbuf,
                   qfile.main_buf,
                   qfile.header_len);

    char *info_s;

    asprintf(&info_s,
             "##INFO=<ID=N_CASE,Number=1,Type=Integer,"
             "Description=\"Number of cases\">\n");
    append_out_buf(&outbuf, info_s, strlen(info_s));

    asprintf(&info_s,
             "##INFO=<ID=N_CTRL,Number=1,Type=Integer,"
             "Description=\"Number of controls\">\n");
    append_out_buf(&outbuf, info_s, strlen(info_s));

    asprintf(&info_s,
             "##INFO=<ID=O_CASE,Number=1,Type=Integer,"
             "Description=\"Number of variants observed in cases\">\n");
    append_out_buf(&outbuf, info_s, strlen(info_s));

    asprintf(&info_s,
             "##INFO=<ID=O_CTRL,Number=1,Type=Integer,"
             "Description=\"Number of variants observed in controls\">\n");
    append_out_buf(&outbuf, info_s, strlen(info_s));

    /*
    asprintf(&info_s,
             "##INFO=<ID=P_CASE_CTRL,Number=%u,Type=Integer,"
             "Description=\"Number of variants in permuted cases "
             "and controls where i is the case and i+1 is the control\">\n",
             N*2);
    append_out_buf(&outbuf, info_s, strlen(info_s));
    */


    char last_header_line[]="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    append_out_buf(&outbuf, last_header_line, strlen(last_header_line));

    for (i=0; i < num_variants; ++i) {

        append_out_buf(&outbuf,
                       qfile.lines[i],
                       qfile.line_lens[i]-1);

        /*
        P_csv_i = 0;
        // The first two values for R are the observations, the rest are 
        // permutation case/ctrl pairs
        for (j = 2; j < N*2 + 2; j+=2) {
            if (j == 2)
                l = sprintf(P_csv + P_csv_i,
                            "%d,%d",
                            R[i*(N*2 + 2) + j],
                            R[i*(N*2 + 2) + j+1]);
            else
                l = sprintf(P_csv + P_csv_i,
                            ",%d,%d",
                            R[i*(N*2 + 2) + j],
                            R[i*(N*2 + 2) + j+1]);
            P_csv_i += l;
        }
        */

        asprintf(&info_s,
                ";N_CASE=%u;"
                "N_CTRL=%u;"
                "O_CASE=%u;"
                "O_CTRL=%u\n",
                //"P_CASE_CTRL=%s\n",
                n_cases,
                n_ctrls,
                R[i*(N*2+2) + 0],
                R[i*(N*2+2) + 1]);
                //P_csv);
        append_out_buf(&outbuf, info_s, strlen(info_s));
    }
    free(info_s);
    quick_file_delete(&qfile);
    free_out_buf(&outbuf);
}
//}}}

//{{{int pop_help(char *func)
int pop_help(char *func)
{
    fprintf(stderr,
"%s %s\n"
"usage:   gqt %s -i <gqt file> \\\n"
"                   -d <ped database file> \\\n"
"                   -f <label db field name> (requried for pca-shared)\\\n"
"                   -l <label output file> (requried for pca-shared)\\\n"
"                   -p <population query 1> \\\n"
"                   -p <population query 2> \\\n"
"                   ... \\\n"
"                   -p <population query N> \n"
"\n"
"Each population query defines one subpopulation.\n"
"For example, to find compare the GBR and YRI subpopulations:\n"
"\t-p \"Population = 'GBR'\"\n"
"\t-p \"Population = 'YRI'\"\n"
"\n"
"Population queries are based on the PED file that is associated with the\n"
"genotypes, and any column in that PED file can be part of the query.  For\n"
"example, a PED file that includes the \"Paternal_ID\" and \"Gender\" fields\n"
"(where male = 1 and female = 2) could be queried by:\n"
"\n"
"\t-p \"Paternal_ID = 'NA12878' AND Gender = 2\"\n"
"\n"
"NOTE: gst and fst assume that variants are biallelic.  If your data\n"
"contains multiallelic sites, we recommend decomposing your VCF \n"
"(see A. Tan, Bioinformatics 2015) prior to indexing.\n",
            PROGRAM_NAME, VERSION, func);
    return EX_USAGE;
}
//}}}

