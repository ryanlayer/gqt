#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

int sum_help();

int sum_plt(char *in,
              uint32_t *R,
              uint32_t num_records,
              uint32_t range_start,
              uint32_t range_end,
              uint32_t range_exclude,
              int x_is_set,
              int time,
              int quiet,
              char *bim);
int sum_ubin(char *in,
               uint32_t *R,
               uint32_t num_records,
               uint32_t range_start,
               uint32_t range_end,
               uint32_t range_exclude,
               int x_is_set,
               int time,
               int quiet,
               char *bim);
int sum_wah(char *in,
              uint32_t *R,
              uint32_t num_records,
              uint32_t range_start,
              uint32_t range_end,
              uint32_t range_exclude,
              int x_is_set,
              int time,
              int quiet,
              char *bim);
int sum_wahbm(char *in,
                uint32_t *R,
                uint32_t num_records,
                uint32_t range_start,
                uint32_t range_end,
                uint32_t range_exclude,
                int x_is_set,
                int time,
                int quiet,
                char *bim);

int sum_in_place_wahbm(char *in,
                         uint32_t *R,
                         uint32_t num_records,
                         uint32_t range_start,
                         uint32_t range_end,
                         uint32_t range_exclude,
                         int x_is_set,
                         int a_is_set,
                         int time,
                         int quiet,
                         char *bim);

int sum_compressed_in_place_wahbm(char *in,
                                    uint32_t *R,
                                    uint32_t num_records,
                                    uint32_t range_start,
                                    uint32_t range_end,
                                    uint32_t range_exclude,
                                    int x_is_set,
                                    int time,
                                    int quiet,
                                    char *bim);


void print_sum_result(uint32_t *R,
                        uint32_t num_fields,
                        char *bim);


int sum(int argc, char **argv)
{
    if (argc < 2) return sum_help();

    int avx_is_on = 0;

#ifdef __AVX2__
    avx_is_on = 1;
#endif

    int c;
    char *in, *out, *record_ids, *op, *bim, *query, *ped_db_file;
    uint32_t num_records;
    int i_is_set = 0,
        a_is_set = 0,
        b_is_set = 0,
        r_is_set = 0,
        n_is_set = 0,
        Q_is_set = 0,
        t_is_set = 0,
        u_is_set = 0,
        l_is_set = 0,
        e_is_set = 0,
        q_is_set = 0,
        d_is_set = 0,
        x_is_set = 0;

    uint32_t range_start = 0, range_end = 3, range_exclude = 0;

    while ((c = getopt (argc, argv, "ahi:r:n:b:Qtu:l:e:x:q:d:")) != -1) {
        switch (c) {
        case 'a':
            a_is_set = 1;
            break;
        case 'd':
        	d_is_set = 1;
        	ped_db_file = optarg;
        	break;
        case 'q':
        	q_is_set = 1;
        	query = optarg;
        	break;
        case 'b':
        	b_is_set = 1;
        	bim = optarg;
        	break;
        case 'u':
            u_is_set = 1;
            range_end = atoi(optarg);
            break;
        case 'l':
            l_is_set = 1;
            range_start = atoi(optarg);
            break;
        case 'x':
            x_is_set = 1;
            range_exclude = atoi(optarg);
            break;
        case 't':
            t_is_set = 1;
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
        case 'Q':
            Q_is_set = 1;
            break;
        case 'h':
            return sum_help();
        case '?':
            if (    (optopt == 'i') ||
                    (optopt == 'r') ||
                    (optopt == 'n') ||
                    (optopt == 'g') ||
                    (optopt == 'l') ||
                    (optopt == 'e') ||
                    (optopt == 'x'))
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

    if ( ((q_is_set == 1) && (d_is_set == 0)) ||
         ((q_is_set == 0) && (d_is_set == 1)) ) {
        printf("Both query and PED database name must be set\n");
        return sum_help();
    }

    if ((q_is_set == 0) && (d_is_set == 0)) {
        if (n_is_set == 0) {
            printf("Number of records is not set\n");
            return sum_help();
        }

        if (r_is_set == 0) {
            printf("Record IDs are not set\n");
            return sum_help();
        }
    }

    uint32_t *R;

    if (q_is_set == 1) {
        num_records = resolve_ind_query(&R, query, ped_db_file);
        //fprintf(stderr, "num_records:%u\n", num_records);
    } else {
        R = (uint32_t *) malloc(num_records * sizeof(uint32_t));
        parse_cmd_line_int_csv(R, num_records, record_ids);
    }

    if ((a_is_set == 1) && (avx_is_on == 0 )) {
        printf("AVX support not included at compile time\n");
        return sum_help();
    }

    if (strcmp(type, "plt") == 0)
        return sum_plt(in,
                         R,
                         num_records,
                         range_start,
                         range_end,
                         range_exclude,
                         x_is_set,
                         t_is_set,
                         Q_is_set,
                         bim);

    else if (strcmp(type, "ubin") == 0)
        return sum_ubin(in,
                          R,
                          num_records,
                          range_start,
                          range_end,
                          range_exclude,
                          x_is_set,
                          t_is_set,
                          Q_is_set,
                          bim);

    else if (strcmp(type, "wah") == 0)
        return sum_wah(in,
                         R,
                         num_records,
                         range_start,
                         range_end,
                         range_exclude,
                         x_is_set,
                         t_is_set,
                         Q_is_set,
                         bim);

    else if (strcmp(type, "wahbm") == 0)
        return sum_wahbm(in,
                           R,
                           num_records,
                           range_start,
                           range_end,
                           range_exclude,
                           x_is_set,
                           t_is_set,
                           Q_is_set,
                           bim);

    else if (strcmp(type, "ipwahbm") == 0)
        return sum_in_place_wahbm(in,
                                    R,
                                    num_records,
                                    range_start,
                                    range_end,
                                    range_exclude,
                                    x_is_set,
                                    a_is_set,
                                    t_is_set,
                                    Q_is_set,
                                    bim);

    else if (strcmp(type, "cipwahbm") == 0)
        return sum_compressed_in_place_wahbm(in,
                                             R,
                                             num_records,
                                             range_start,
                                             range_end,
                                             range_exclude,
                                             x_is_set,
                                             t_is_set,
                                             Q_is_set,
                                             bim);


    return 1;
}

int sum_help()
{
    printf("usage:   gqt sum <type> -i <input file>\n"
           "                        -n <number of records>\n"
           "                        -r <record ids CSV>\n"
           "                        -l <inclusive lower bound>\n"
           "                        -u <inclusive upper bound>\n"
           "                        -b <bim file>\n"
           "                        -x <exlcude>\n"
           "                        -q <query>\n"
           "                        -d <ped database>\n"
           "\ttypes:\n"
           "\t\tplt       Plain text \n"
           "\t\tubin      Uncompressed binary\n"
           "\t\twah       WAH \n"
           "\t\twahbm     WAH bitmap\n"
           "\t\tipwahbm   in-place WAH bitmap\n"
           "\t\tcipwahbm  compressed in-place WAH bitmap\n"
          );

    return 0;
}

int sum_plt(char *in,
              uint32_t *R,
              uint32_t num_records,
              uint32_t range_start,
              uint32_t range_end,
              uint32_t range_exclude,
              int x_is_set,
              int time,
              int quiet,
              char *bim)
{
#if 0
    start();
    struct plt_file pf = init_plt_file(in);
    uint32_t *pf_R;

    uint32_t len_pf_R;
    
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
        print_sum_result(pf_R, pf.num_fields, bim);


    free(pf_R);
    fclose(pf.file);
#endif
    return 0;
}

int sum_ubin(char *in,
               uint32_t *R,
               uint32_t num_records,
               uint32_t range_start,
               uint32_t range_end,
               uint32_t range_exclude,
               int x_is_set,
               int time,
               int quiet,
               char *bim)

{
#if 0
    start();
    struct ubin_file uf = init_ubin_file(in);
    uint32_t *uf_R;
    uint32_t len_uf_R;

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
        print_sum_result(uf_R, uf.num_fields, bim);

    free(uf_R);
    fclose(uf.file);
#endif
    return 0;
}

int sum_wah(char *in,
              uint32_t *R,
              uint32_t num_records,
              uint32_t range_start,
              uint32_t range_end,
              uint32_t range_exclude,
              int x_is_set,
              int time,
              int quiet,
              char *bim)

{
    return 0;
}

int sum_in_place_wahbm(char *in,
                         uint32_t *R,
                         uint32_t num_records,
                         uint32_t range_start,
                         uint32_t range_end,
                         uint32_t range_exclude,
                         int x_is_set,
                         int a_is_set,
                         int time,
                         int quiet,
                         char *bim)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R;

    if (time != 0 )
        start();

    if (a_is_set ==0)
        len_wf_R = sum_range_records_in_place_wahbm(wf,
                                                R,
                                                num_records,
                                                range_start,
                                                range_end + 1,
                                                &wf_R);
#ifdef __AVX2__
    else
        len_wf_R = avx_sum_range_records_in_place_wahbm(wf,
                                                R,
                                                num_records,
                                                range_start,
                                                range_end + 1,
                                                &wf_R);

#endif

    if (time != 0 ) {
        stop();
        fprintf(stderr,"%lu\n", report());
    }

    if (quiet == 0)
        print_sum_result(wf_R, wf.num_fields, bim);

    free(wf_R);
    fclose(wf.file);

    return 0;

}

int sum_compressed_in_place_wahbm(char *in,
                                    uint32_t *R,
                                    uint32_t num_records,
                                    uint32_t range_start,
                                    uint32_t range_end,
                                    uint32_t range_exclude,
                                    int x_is_set,
                                    int time,
                                    int quiet,
                                    char *bim)

{
#if 0
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R;

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
        print_sum_result(wf_R, wf.num_fields, bim);

    free(wf_R);
    fclose(wf.file);
#endif
    return 0;



}

int sum_wahbm(char *in,
                uint32_t *R,
                uint32_t num_records,
                uint32_t range_start,
                uint32_t range_end,
                uint32_t range_exclude,
                int x_is_set,
                int time,
                int quiet,
                char *bim)

{
#if 0
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R;

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
        print_sum_result(wf_R, wf.num_fields, bim);

    free(wf_R);
    fclose(wf.file);
#endif
    return 0;
}

void print_sum_result(uint32_t *R,
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
