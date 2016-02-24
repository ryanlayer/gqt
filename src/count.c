#ifdef __AVX__
#include <immintrin.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#include "wahbm.h"
#include "wahbm_in_place.h"
#include "wahbm_compressed_in_place.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

int count_help();

int count_plt(char *in,
              uint32_t query_value,
              char *op,
              uint32_t *R,
              uint32_t num_records,
              int time,
              int quiet,
              char *bim);

int count_ubin(char *in,
               uint32_t query_value,
               char *op,
               uint32_t *R,
               uint32_t num_records,
               int time,
               int quiet,
               char *bim);
int count_wah(char *in,
              uint32_t query_value,
              char *op,
              uint32_t *R,
              uint32_t num_records,
              int time,
              int quiet,
              char *bim);
int count_wahbm(char *in,
                uint32_t query_value,
                char *op,
                uint32_t *R,
                uint32_t num_records,
                int time,
                int quiet,
                char *bim);

int count_in_place_wahbm(char *in,
                         uint32_t query_value,
                         char *op,
                         uint32_t *R,
                         uint32_t num_records,
                         int time,
                         int quiet,
                         int avx,
                         char *bim);

int count_compressed_in_place_wahbm(char *in,
                                    uint32_t query_value,
                                    char *op,
                                    uint32_t *R,
                                    uint32_t num_records,
                                    int time,
                                    int quiet,
                                    char *bim);


void print_count_result(uint32_t *R,
                        uint32_t num_fields,
                        char *bim);


int count(int argc, char **argv)
{
    if (argc < 2) return count_help();

    int c;
    char *in=NULL, *out=NULL, *record_ids=NULL, *op=NULL, *bim=NULL;
    uint32_t query_value = 0, num_records = 0;
    int i_is_set = 0,
        a_is_set = 0,
        b_is_set = 0,
        o_is_set = 0,
        r_is_set = 0,
        n_is_set = 0,
        Q_is_set = 0,
        t_is_set = 0,
        q_is_set = 0;

    while ((c = getopt (argc, argv, "hi:o:q:r:n:b:Qta")) != -1) {
        switch (c) {
        case 'a':
            a_is_set = 1;
            break;
        case 'b':
        	b_is_set = 1;
        	bim = optarg;
        	break;
        case 't':
            t_is_set = 1;
            break;
        case 'o':
            o_is_set = 1;
            op = optarg;
            break;
        case 'r':
            r_is_set = 1;
            record_ids= optarg;
            break;
        case 'n':
            n_is_set = 1;
            num_records = atoi(optarg);
            break;
        case 'i':
            i_is_set = 1;
            in = optarg;
            break;
        case 'q':
            q_is_set = 1;
            query_value = atoi(optarg);
            break;
        case 'Q':
            Q_is_set = 1;
            break;
        case 'h':
            return count_help();
        case '?':
            if ( (optopt == 'i') ||
                    (optopt == 'o') ||
                    (optopt == 'r') ||
                    (optopt == 'n') ||
                    (optopt == 'q') )
                fprintf (stderr, "Option -%c requires an argument.\n",
                         optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        default:
            return count_help();
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return count_help();
    }

    if (q_is_set == 0) {
        printf("Query value is not set\n");
        return count_help();
    }

    if (n_is_set == 0) {
        printf("Number of records is not set\n");
        return count_help();
    }

    if (r_is_set == 0) {
        printf("Record IDs are not set\n");
        return count_help();
    }

    if (o_is_set == 0) {
        printf("Opperation is not set\n");
        return count_help();
    }


    if ( !((strcmp(op, "gt") == 0) ||
           (strcmp(op, "lt") == 0) ||
           (strcmp(op, "eq") == 0) ||
           (strcmp(op, "ne") == 0) ||
           (strcmp(op, "le") == 0) ||
           (strcmp(op, "ge") == 0)) ) {
        printf("Unknown opperation\n");
        return count_help();
    }


    uint32_t R[num_records];
    parse_cmd_line_int_csv(R, num_records, record_ids);

    if (strcmp(type, "plt") == 0)
        return count_plt(in,
                         query_value,
                         op,
                         R,
                         num_records,
                         t_is_set,
                         Q_is_set,
                         bim);

    else if (strcmp(type, "ubin") == 0)
        return count_ubin(in,
                          query_value,
                          op,
                          R,
                          num_records,
                          t_is_set,
                          Q_is_set,
                          bim);

    else if (strcmp(type, "wah") == 0)
        return count_wah(in,
                         query_value,
                         op,
                         R,
                         num_records,
                         t_is_set,
                         Q_is_set,
                         bim);

    else if (strcmp(type, "wahbm") == 0)
        return count_wahbm(in,
                           query_value,
                           op,
                           R,
                           num_records,
                           t_is_set,
                           Q_is_set,
                           bim);

    else if (strcmp(type, "ipwahbm") == 0)
        return count_in_place_wahbm(in,
                                    query_value,
                                    op,
                                    R,
                                    num_records,
                                    t_is_set,
                                    Q_is_set,
                                    a_is_set,
                                    bim);

    else if (strcmp(type, "cipwahbm") == 0)
        return count_compressed_in_place_wahbm(in,
                                               query_value,
                                               op,
                                               R,
                                               num_records,
                                               t_is_set,
                                               Q_is_set,
                                               bim);


    return 1;
}

int count_help()
{
    printf("usage:   gqt count <type> -o <opperation> -i <input file> "
           "-q <query value> -n <number of records> -r <record ids> -b<bim file>\n"
           "op:\n"
           "        gt         Greater than\n"
           "        lt         Less than\n"
           "        eq         Equal\n"
           "        ne         Not equal\n"
           "        le         Less than or equal\n"
           "        ge         Greater than or equal\n"
           "type:"
           "         plt       Plain text \n"
           "         ubin      Uncompressed binary\n"
           "         wah       WAH \n"
           "         wahbm     WAH bitmap\n"
           "         ipwahbm   in-place WAH bitmap\n"
           "         cipwahbm  compressed in-place WAH bitmap\n"
          );

    return 0;
}

int count_plt(char *in,
              uint32_t query_value,
              char *op,
              uint32_t *R,
              uint32_t num_records,
              int time,
              int quiet,
              char *bim)
{
    start();
    struct plt_file pf = init_plt_file(in);
    uint32_t *pf_R;

    uint32_t len_pf_R;
    
    if (strcmp(op,"gt") == 0)
        len_pf_R = gt_count_records_plt(pf,
                                        R,
                                        num_records,
                                        query_value,
                                        &pf_R);
    else 
        return count_help();

    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_count_result(pf_R, pf.num_fields, bim);


    free(pf_R);
    fclose(pf.file);

    return 0;
}

int count_ubin(char *in,
               uint32_t query_value,
               char *op,
               uint32_t *R,
               uint32_t num_records,
               int time,
               int quiet,
               char *bim)

{
    start();
#if 0
    struct ubin_file uf = init_ubin_file(in);
    uint32_t *uf_R;
    uint32_t len_uf_R;

    if (strcmp(op,"gt") == 0)
        len_uf_R = gt_count_records_ubin(uf,
                                         R,
                                         num_records,
                                         query_value,
                                         &uf_R);
    else 
        return count_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_count_result(uf_R, uf.num_fields);

    free(uf_R);
    fclose(uf.file);
#endif
    return 0;
}

int count_wah(char *in,
              uint32_t query_value,
              char *op,
              uint32_t *R,
              uint32_t num_records,
              int time,
              int quiet,
              char *bim)

{
    return 0;
}

int count_in_place_wahbm(char *in,
                         uint32_t query_value,
                         char *op,
                         uint32_t *R,
                         uint32_t num_records,
                         int time,
                         int quiet,
                         int avx,
                         char *bim)

{
    if (time != 0 )
        start();

    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);
    __attribute__((aligned(64)))uint32_t *wf_R = NULL; 
    //__declspec(align(64)) uint32_t *wf_R;
    uint32_t len_wf_R;

    if (strcmp(op,"gt") == 0) {
        if (avx == 0)
            len_wf_R = gt_count_records_in_place_wahbm(wf,
                                                       R,
                                                       num_records,
                                                       query_value,
                                                       &wf_R);
#ifdef __AVX2__
        else
            len_wf_R = avx_gt_count_records_in_place_wahbm(wf,
                                                       R,
                                                       num_records,
                                                       query_value,
                                                       &wf_R);
#endif

    } else 
        return count_help();

    if (time != 0 ) {
        stop();
        fprintf(stderr,"%lu\n", report());
    }

    if (quiet == 0)
        //print_count_result(wf_R, wf.num_fields, bim);
        print_count_result(wf_R, wf->gqt_header->num_variants, bim);

    free(wf_R);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);
    return 0;
}

int count_compressed_in_place_wahbm(char *in,
                                    uint32_t query_value,
                                    char *op,
                                    uint32_t *R,
                                    uint32_t num_records,
                                    int time,
                                    int quiet,
                                    char *bim)

{
    start();
    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R;

    if (strcmp(op,"gt") == 0)
        len_wf_R = gt_count_records_compressed_in_place_wahbm(wf,
                                                              R,
                                                              num_records,
                                                              query_value,
                                                              &wf_R);
    else 
        return count_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_count_result(wf_R, wf->gqt_header->num_variants, bim);

    free(wf_R);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}

int count_wahbm(char *in,
                uint32_t query_value,
                char *op,
                uint32_t *R,
                uint32_t num_records,
                int time,
                int quiet,
                char *bim)

{
    start();
    //struct wah_file wf = init_wahbm_file(in);
    struct wahbm_file *wf = open_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R;

    if (strcmp(op,"gt") == 0)
        len_wf_R = gt_count_records_wahbm(wf,
                                          R,
                                          num_records,
                                          query_value,
                                          &wf_R);
    else 
        return count_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        //print_count_result(wf_R, wf.num_fields, bim);
        print_count_result(wf_R, wf->gqt_header->num_variants, bim);

    free(wf_R);
    //destroy_wah_file(&wf);
    destroy_wahbm_file(wf);

    return 0;
}

void print_count_result(uint32_t *R,
                        uint32_t num_fields,
                        char *bim)
{

    struct quick_file_info qfile;
    struct output_buffer out_buf;
    size_t i=0;

    if (bim != NULL) {
    	quick_file_init(bim, &qfile);
    }

    init_out_buf(&out_buf, NULL);


    for(; i < num_fields; ++i) {
        if (bim != NULL) {
	    append_out_buf(&out_buf, qfile.lines[i], qfile.line_lens[i]);
        }
        append_out_buf(&out_buf, "\t", 1);
        append_integer_to_out_buf(&out_buf,  R[i]);
        append_out_buf(&out_buf, "\n", 1);
    }
    append_out_buf(&out_buf, "\n", 1);
    quick_file_delete(&qfile);
    free_out_buf(&out_buf);
}
