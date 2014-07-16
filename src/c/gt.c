#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"
#include "timer.h"

int gt_help();

int gt_plt(char *in,
           unsigned int query_value,
           unsigned int *R,
           unsigned int num_records,
           int time,
           int quiet);

int gt_ubin(char *in,
            unsigned int query_value,
            unsigned int *R,
            unsigned int num_records,
            int time,
            int quiet);

int gt_wah(char *in,
           unsigned int query_value,
           unsigned int *R,
           unsigned int num_records,
           int time,
           int quiet);

int gt_wahbm(char *in,
             unsigned int query_value,
             unsigned int *R,
             unsigned int num_records,
             int time,
             int quiet);

int gt_in_place_wahbm(char *in,
                      unsigned int query_value,
                      unsigned int *R,
                      unsigned int num_records,
                      int time,
                      int quiet);

int gt_compressed_in_place_wahbm(char *in,
                                 unsigned int query_value,
                                 unsigned int *R,
                                 unsigned int num_records,
                                 int time,
                                 int quiet);


void print_result(unsigned int len,
                  unsigned int *R,
                  unsigned int num_fields);



int gt(int argc, char **argv)
{
    if (argc < 2) return gt_help();

    int c;
    char *in, *out, *record_ids;
    unsigned int query_value, num_records;
    int i_is_set = 0, 
        r_is_set = 0, 
        n_is_set = 0, 
        Q_is_set = 0, 
        t_is_set = 0, 
        q_is_set = 0; 

    while ((c = getopt (argc, argv, "hi:q:r:n:Qt")) != -1) {
        switch (c) {
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
            case 't':
                t_is_set = 1;
                break;
            case 'h':
                gt_help();
                return 1;
            case '?':
                if ( (optopt == 'i') || (optopt == 'q') )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                gt_help();
                return 1;
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return gt_help();
    } 

    if (q_is_set == 0) {
        printf("Query value is not set\n");
        return gt_help();
    } 

    if (n_is_set == 0) {
        printf("Number of records is not set\n");
        return gt_help();
    } 

    if (r_is_set == 0) {
        printf("Record IDs are not set\n");
        return gt_help();
    } 

    unsigned int R[num_records];
    parse_cmd_line_int_csv(R, num_records, record_ids);

    if (strcmp(type, "plt") == 0)
        return gt_plt(in,
                      query_value,
                      R,
                      num_records,
                      t_is_set,
                      Q_is_set);

    else if (strcmp(type, "ubin") == 0)
        return gt_ubin(in,
                       query_value,
                       R,
                       num_records,
                      t_is_set,
                       Q_is_set);

    else if (strcmp(type, "wah") == 0) 
        return gt_wah(in,
                      query_value,
                      R,
                      num_records,
                      t_is_set,
                      Q_is_set);

    else if (strcmp(type, "wahbm") == 0)
        return gt_wahbm(in,
                        query_value,
                        R,
                        num_records,
                        t_is_set,
                        Q_is_set);

    else if (strcmp(type, "ipwahbm") == 0)
        return gt_in_place_wahbm(in,
                                 query_value,
                                 R,
                                 num_records,
                                 t_is_set,
                                 Q_is_set);

    else if (strcmp(type, "cipwahbm") == 0)
        return gt_compressed_in_place_wahbm(in,
                                            query_value,
                                            R,
                                            num_records,
                                            t_is_set,
                                            Q_is_set);




    return 1;
}

int gt_help()
{
    printf("usage:   gtq gt <type> -i <input file> -q <query value> "
                "-n <number of records> -r <record ids>\n"
           "         plt       Plain text \n"
           "         ubin      Uncompressed binary\n"
           "         wah       WAH \n"
           "         wahbm     WAH bitmap\n"
           "         ipwahbm   in-place WAH bitmap\n"
           "         cipwahbm  compressed in-place WAH bitmap\n"
    );

    return 0;
}

int gt_plt(char *in,
           unsigned int query_value,
           unsigned int *R,
           unsigned int num_records,
           int time,
           int quiet)
{
    start();
    struct plt_file pf = init_plt_file(in);
    unsigned int *pf_R;
    unsigned int len_pf_R = gt_records_plt(pf,
                                           R,
                                           num_records,
                                           query_value,
                                           &pf_R);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_pf_R, pf_R, pf.num_fields);


    free(pf_R);
    fclose(pf.file);

    return 0;
}

int gt_ubin(char *in,
            unsigned int query_value,
            unsigned int *R,
            unsigned int num_records,
            int time,
            int quiet)

{
    start();
    struct ubin_file uf = init_ubin_file(in);
    unsigned int *uf_R;
    unsigned int len_uf_R = gt_records_ubin(uf,
                                           R,
                                           num_records,
                                           query_value,
                                           &uf_R);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_uf_R, uf_R, uf.num_fields);

    free(uf_R);
    fclose(uf.file);

    return 0;
}
int gt_wah(char *in,
           unsigned int query_value,
           unsigned int *R,
           unsigned int num_records,
           int time,
           int quiet)

{
    return 0;
}

int gt_in_place_wahbm(char *in,
                      unsigned int query_value,
                      unsigned int *R,
                      unsigned int num_records,
                      int time,
                      int quiet)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_in_place_wahbm(wf,
                                                      R,
                                                      num_records,
                                                      query_value,
                                                      &wf_R);
    unsigned int *ints;
    unsigned int len_ints = wah_to_ints(wf_R,len_wf_R,&ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

int gt_compressed_in_place_wahbm(char *in,
                                 unsigned int query_value,
                                 unsigned int *R,
                                 unsigned int num_records,
                                 int time,
                                 int quiet)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_compressed_in_place_wahbm(wf,
                                                                 R,
                                                                 num_records,
                                                                 query_value,
                                                                 &wf_R);
    unsigned int *ints;
    unsigned int len_ints = compressed_in_place_wah_to_ints(wf_R,
                                                            len_wf_R,
                                                            &ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

int gt_wahbm(char *in,
             unsigned int query_value,
             unsigned int *R,
             unsigned int num_records,
             int time,
             int quiet)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_wahbm(wf,
                                             R,
                                             num_records,
                                             query_value,
                                             &wf_R);

    unsigned int *ints;
    unsigned int len_ints = wah_to_ints(wf_R,len_wf_R,&ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

void print_result(unsigned int len,
                  unsigned int *R,
                  unsigned int num_fields)
{
    unsigned int i,j, bit_i = 0;
    for(i = 0; i < len; ++i) {
        if (i!= 0)
            printf(" ");

        int *r = unpack_1_bit_ints(R[i]);

        for(j = 0; j < 32; ++j) {
            if (j!= 0)
                printf(" ");
            printf("%d", r[j]);

            bit_i += 1;
            if (bit_i == num_fields)
                break;
        }
        if (bit_i == num_fields)
            break;
        free(r);
    }
    printf("\n");
}
