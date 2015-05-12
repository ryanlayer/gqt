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

int pop_help(char *op);


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

void print_pop_result(char *op,
                      double *R,
                      uint32_t num_variants,
                      char *bim);

#if 0
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
#endif

//{{{ int pop(char *op, int argc, char **argv)
int pop(char *op, int argc, char **argv)
{
    if (argc < 2) return pop_help(op);

    int c;
    char *wahbm_file_name=NULL,
         *id_query=NULL,
         *db_file_name=NULL,
         *bim_file_name=NULL,
         *vid_file_name=NULL;
    int i_is_set = 0,
        id_q_count = 0,
        d_is_set = 0,
        v_is_set = 0,
        b_is_set = 0;

    char *id_query_list[100];

    //{{{ parse cmd line opts
    while ((c = getopt (argc, argv, "hi:p:d:b:v:")) != -1) {
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
                 (optopt == 'p') ||
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
            printf("Auto detect failure: VID file %s not found\n",
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
            printf("Auto detect failure: PED DB file %s not found\n",
                   db_file_name);
            return pop_help(op);
        }
    }

    if (i_is_set == 0) {
        printf("GQT file is not set\n");
        return pop_help(op);
    }

    if (v_is_set == 0) {
        printf("VID file is not set\n");
        return pop_help(op);
    }

    if (b_is_set == 0) {
        printf("BIM file is not set\n");
        return pop_help(op);
    }

    if (d_is_set == 0) {
        printf("PED database file is not set\n");
        return pop_help(op);
    }
    //}}}

    // open WAH/GQT file
    struct wah_file wf = init_wahbm_file(wahbm_file_name);

    // open VID file
    FILE *vid_f = fopen(vid_file_name, "rb");
    if (vid_f == NULL) {
        fprintf(stderr, "Could not read VIDE file: %s\n", vid_file_name);
        return 1;
    }
    uint32_t *vids = (uint32_t *) malloc(wf.num_fields*sizeof(uint32_t));
    int r = fread(vids, sizeof(uint32_t), wf.num_fields, vid_f);
    fclose(vid_f);

    uint32_t num_ints = (wf.num_fields + 32 - 1)/ 32;
    uint32_t len_ints;

    double *R;
    uint32_t len_R;

    if (strcmp("fst",op) == 0) {
        len_R = fst(&wf, id_query_list, id_q_count, vids, db_file_name, &R);
    } else if (strcmp("gst",op) == 0) {
        len_R = gst(&wf, id_query_list, id_q_count, vids, db_file_name, &R);
    } else {
        fprintf(stderr, "ERROR: Unknown opperation %s", op);
        pop_help("fst|gst");
    }

    print_pop_result(op, R, wf.num_fields, bim_file_name);

    fclose(wf.file);
    return 0;
}

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

//{{{int pop_help(char *func)
int pop_help(char *func)
{
    printf(
"usage:   gqt %s -i <gqt file> \\\n"
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
"\t-p \"Paternal_ID = 'NA12878' AND Gender = 2\"\n", func);
    return 1;
}
//}}}
