#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"
#include "timer.h"

int sum_help();

int sum_plt(char *in,
              unsigned int query_value,
              char *op,
              unsigned int *R,
              unsigned int num_records,
              int time,
              int quiet);
int sum_ubin(char *in,
               unsigned int query_value,
               char *op,
               unsigned int *R,
               unsigned int num_records,
               int time,
               int quiet);
int sum_wah(char *in,
              unsigned int query_value,
              char *op,
              unsigned int *R,
              unsigned int num_records,
              int time,
              int quiet);
int sum_wahbm(char *in,
                unsigned int query_value,
                char *op,
                unsigned int *R,
                unsigned int num_records,
                int time,
                int quiet);

int sum_in_place_wahbm(char *in,
                         unsigned int query_value,
                         char *op,
                         unsigned int *R,
                         unsigned int num_records,
                         int time,
                         int quiet);

int sum_compressed_in_place_wahbm(char *in,
                                    unsigned int query_value,
                                    char *op,
                                    unsigned int *R,
                                    unsigned int num_records,
                                    int time,
                                    int quiet);


void print_sum_result(unsigned int *R,
                        unsigned int num_fields);


int sum(int argc, char **argv)
{
    if (argc < 2) return sum_help();

    int c;
    char *in, *out, *record_ids, *op;
    unsigned int query_value, num_records;
    int i_is_set = 0,
        o_is_set = 0,
        r_is_set = 0,
        n_is_set = 0,
        Q_is_set = 0,
        t_is_set = 0,
        q_is_set = 0;

    while ((c = getopt (argc, argv, "hi:o:q:r:n:Qt")) != -1) {
        switch (c) {
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
            return sum_help();
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
            return sum_help();
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return sum_help();
    }

    if (q_is_set == 0) {
        printf("Query value is not set\n");
        return sum_help();
    }

    if (n_is_set == 0) {
        printf("Number of records is not set\n");
        return sum_help();
    }

    if (r_is_set == 0) {
        printf("Record IDs are not set\n");
        return sum_help();
    }

    if (o_is_set == 0) {
        printf("Opperation is not set\n");
        return sum_help();
    }


    if ( !((strcmp(op, "gt") == 0) ||
           (strcmp(op, "lt") == 0) ||
           (strcmp(op, "eq") == 0) ||
           (strcmp(op, "ne") == 0) ||
           (strcmp(op, "le") == 0) ||
           (strcmp(op, "ge") == 0)) ) {
        printf("Unknown opperation\n");
        return sum_help();
    }


    unsigned int R[num_records];
    parse_cmd_line_int_csv(R, num_records, record_ids);

    if (strcmp(type, "plt") == 0)
        return sum_plt(in,
                         query_value,
                         op,
                         R,
                         num_records,
                         t_is_set,
                         Q_is_set);

    else if (strcmp(type, "ubin") == 0)
        return sum_ubin(in,
                          query_value,
                          op,
                          R,
                          num_records,
                          t_is_set,
                          Q_is_set);

    else if (strcmp(type, "wah") == 0)
        return sum_wah(in,
                         query_value,
                         op,
                         R,
                         num_records,
                         t_is_set,
                         Q_is_set);

    else if (strcmp(type, "wahbm") == 0)
        return sum_wahbm(in,
                           query_value,
                           op,
                           R,
                           num_records,
                           t_is_set,
                           Q_is_set);

    else if (strcmp(type, "ipwahbm") == 0)
        return sum_in_place_wahbm(in,
                                    query_value,
                                    op,
                                    R,
                                    num_records,
                                    t_is_set,
                                    Q_is_set);

    else if (strcmp(type, "cipwahbm") == 0)
        return sum_compressed_in_place_wahbm(in,
                                               query_value,
                                               op,
                                               R,
                                               num_records,
                                               t_is_set,
                                               Q_is_set);


    return 1;
}

int sum_help()
{
    printf("usage:   gtq sum <type> -o <opperation> -i <input file> "
           "-q <query value> -n <number of records> -r <record ids>\n"
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

int sum_plt(char *in,
              unsigned int query_value,
              char *op,
              unsigned int *R,
              unsigned int num_records,
              int time,
              int quiet)
{
#if 0
    start();
    struct plt_file pf = init_plt_file(in);
    unsigned int *pf_R;

    unsigned int len_pf_R;
    
    if (strcmp(op,"gt") == 0)
        len_pf_R = gt_sum_records_plt(pf,
                                        R,
                                        num_records,
                                        query_value,
                                        &pf_R);
    else 
        return sum_help();

    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_sum_result(pf_R, pf.num_fields);


    free(pf_R);
    fclose(pf.file);
#endif
    return 0;
}

int sum_ubin(char *in,
               unsigned int query_value,
               char *op,
               unsigned int *R,
               unsigned int num_records,
               int time,
               int quiet)

{
#if 0
    start();
    struct ubin_file uf = init_ubin_file(in);
    unsigned int *uf_R;
    unsigned int len_uf_R;

    if (strcmp(op,"gt") == 0)
        len_uf_R = gt_sum_records_ubin(uf,
                                         R,
                                         num_records,
                                         query_value,
                                         &uf_R);
    else 
        return sum_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_sum_result(uf_R, uf.num_fields);

    free(uf_R);
    fclose(uf.file);
#endif
    return 0;
}

int sum_wah(char *in,
              unsigned int query_value,
              char *op,
              unsigned int *R,
              unsigned int num_records,
              int time,
              int quiet)

{
    return 0;
}

int sum_in_place_wahbm(char *in,
                         unsigned int query_value,
                         char *op,
                         unsigned int *R,
                         unsigned int num_records,
                         int time,
                         int quiet)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R;

    if (strcmp(op,"gt") == 0)
        len_wf_R = gt_sum_records_in_place_wahbm(wf,
                                                   R,
                                                   num_records,
                                                   query_value,
                                                   &wf_R);
    else 
        return sum_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_sum_result(wf_R, wf.num_fields);

    free(wf_R);
    fclose(wf.file);

    return 0;

}

int sum_compressed_in_place_wahbm(char *in,
                                    unsigned int query_value,
                                    char *op,
                                    unsigned int *R,
                                    unsigned int num_records,
                                    int time,
                                    int quiet)

{
#if 0
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R;

    if (strcmp(op,"gt") == 0)
        len_wf_R = gt_sum_records_compressed_in_place_wahbm(wf,
                                                              R,
                                                              num_records,
                                                              query_value,
                                                              &wf_R);
    else 
        return sum_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_sum_result(wf_R, wf.num_fields);

    free(wf_R);
    fclose(wf.file);
#endif
    return 0;



}

int sum_wahbm(char *in,
                unsigned int query_value,
                char *op,
                unsigned int *R,
                unsigned int num_records,
                int time,
                int quiet)

{
#if 0
    start();
    struct wah_file wf = init_wahbm_file(in);
    unsigned int *wf_R;
    unsigned int len_wf_R;

    if (strcmp(op,"gt") == 0)
        len_wf_R = gt_sum_records_wahbm(wf,
                                          R,
                                          num_records,
                                          query_value,
                                          &wf_R);
    else 
        return sum_help();

    stop();

    if (time != 0 )
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_sum_result(wf_R, wf.num_fields);

    free(wf_R);
    fclose(wf.file);
#endif
    return 0;
}

void print_sum_result(unsigned int *R,
                        unsigned int num_fields)
{
    unsigned int i;
    for(i = 0; i < num_fields; ++i) {
        if (i!= 0)
            printf(" ");
        printf("%u", R[i]);
    }
    printf("\n");
}
