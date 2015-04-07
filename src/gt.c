#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"
#include "timer.h"
#include "quick_file.h"
#include "output_buffer.h"

int gt_help();

int gt_plt(char *in,
           uint32_t query_value,
           uint32_t *R,
           uint32_t num_records,
           int time,
           int quiet,
           char *bim);

int gt_plt_fields(char *in,
                  uint32_t query_value,
                  uint32_t *R,
                  uint32_t num_fields,
                  int time,
                  int quiet,
                  char *bim);

int gt_ubin(char *in,
            uint32_t query_value,
            uint32_t *R,
            uint32_t num_records,
            int time,
            int quiet,
            char *bim);

int gt_wah(char *in,
           uint32_t query_value,
           uint32_t *R,
           uint32_t num_records,
           int time,
           int quiet,
           char *bim);

int gt_wahbm(char *in,
             uint32_t query_value,
             uint32_t *R,
             uint32_t num_records,
             int time,
             int quiet,
             char *bim);

int gt_in_place_wahbm(char *in,
                      uint32_t query_value,
                      uint32_t *R,
                      uint32_t num_records,
                      int time,
                      int quiet,
                      char *bim);

int gt_compressed_in_place_wahbm(char *in,
                                 uint32_t query_value,
                                 uint32_t *R,
                                 uint32_t num_records,
                                 int time,
                                 int quiet,
                                 char *bim);


void print_result(uint32_t len,
                  uint32_t *R,
                  uint32_t num_fields,
                  char *bim);



int gt(int argc, char **argv)
{
    if (argc < 2) return gt_help();

    int c;
    char *in, *out, *record_ids, *bim = NULL;
    uint32_t query_value, num_records;
    int i_is_set = 0, 
        r_is_set = 0, 
        f_is_set = 0, 
        n_is_set = 0, 
        Q_is_set = 0, 
        t_is_set = 0, 
        q_is_set = 0,
    	b_is_set = 0;

    while ((c = getopt (argc, argv, "hi:q:r:n:b:fQt")) != -1) {
        switch (c) {
            case 'b':
                b_is_set = 1;
                bim = optarg;
                break;
            case 'r':
                r_is_set = 1;
                record_ids= optarg;
                break;
            case 'n':
                n_is_set = 1;
                num_records = atoi(optarg);
                break;
            case 'f':
                f_is_set = 1;
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

    if (n_is_set == 0) {
        printf("Bim file not set\n");
        return gt_help();
    } 

    uint32_t R[num_records];
    parse_cmd_line_int_csv(R, num_records, record_ids);

    if (strcmp(type, "plt") == 0) {
        if (f_is_set == 1) 
            return gt_plt_fields(in,
                                 query_value,
                                 R,
                                 num_records,
                                 t_is_set,
                                 Q_is_set,
                                 bim);
        else
            return gt_plt(in,
                          query_value,
                          R,
                          num_records,
                          t_is_set,
                          Q_is_set,
                          bim);

    }

    else if (strcmp(type, "ubin") == 0)
        return gt_ubin(in,
                       query_value,
                       R,
                       num_records,
                      t_is_set,
                       Q_is_set,
                       bim);

    else if (strcmp(type, "wah") == 0) 
        return gt_wah(in,
                      query_value,
                      R,
                      num_records,
                      t_is_set,
                      Q_is_set,
                      bim);

    else if (strcmp(type, "wahbm") == 0)
        return gt_wahbm(in,
                        query_value,
                        R,
                        num_records,
                        t_is_set,
                        Q_is_set,
                        bim);

    else if (strcmp(type, "ipwahbm") == 0)
        return gt_in_place_wahbm(in,
                                 query_value,
                                 R,
                                 num_records,
                                 t_is_set,
                                 Q_is_set,
                                 bim);

    else if (strcmp(type, "cipwahbm") == 0)
        return gt_compressed_in_place_wahbm(in,
                                            query_value,
                                            R,
                                            num_records,
                                            t_is_set,
                                            Q_is_set,
                                            bim);




    return 1;
}

int gt_help()
{
    printf("usage:   gqt gt <type> -i <input file> -q <query value> "
                "-n <number of records> -r <record ids> -b <bim file>\n"
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
           uint32_t query_value,
           uint32_t *R,
           uint32_t num_records,
           int time,
           int quiet,
           char *bim)
{
    start();
    struct plt_file pf = init_plt_file(in);
    uint32_t *pf_R;
    uint32_t len_pf_R = gt_records_plt(pf,
                                           R,
                                           num_records,
                                           query_value,
                                           &pf_R);

    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_pf_R, pf_R, pf.num_fields, bim);


    free(pf_R);
    fclose(pf.file);

    return 0;
}

int gt_plt_fields(char *in,
                  uint32_t query_value,
                  uint32_t *R,
                  uint32_t num_records,
                  int time,
                  int quiet,
                  char *bim)
{
    start();
    struct plt_file pf = init_plt_file(in);
    uint32_t *pf_R;
    uint32_t len_pf_R = gt_fields_plt(pf,
                                          R,
                                          num_records,
                                          query_value,
                                          &pf_R);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_pf_R, pf_R, pf.num_records, bim);


    free(pf_R);
    fclose(pf.file);

    return 0;
}

int gt_ubin(char *in,
            uint32_t query_value,
            uint32_t *R,
            uint32_t num_records,
            int time,
            int quiet,
            char *bim)

{
    start();
    struct ubin_file uf = init_ubin_file(in);
    uint32_t *uf_R;
    uint32_t len_uf_R = gt_records_ubin(uf,
                                           R,
                                           num_records,
                                           query_value,
                                           &uf_R);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_uf_R, uf_R, uf.num_fields, bim);

    free(uf_R);
    fclose(uf.file);

    return 0;
}
int gt_wah(char *in,
           uint32_t query_value,
           uint32_t *R,
           uint32_t num_records,
           int time,
           int quiet,
           char *bim)

{
    return 0;
}

int gt_in_place_wahbm(char *in,
                      uint32_t query_value,
                      uint32_t *R,
                      uint32_t num_records,
                      int time,
                      int quiet,
                      char *bim)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R = gt_records_in_place_wahbm(wf,
                                                      R,
                                                      num_records,
                                                      query_value,
                                                      &wf_R);
    uint32_t *ints;
    uint32_t len_ints = wah_to_ints(wf_R,len_wf_R,&ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields, bim);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

int gt_compressed_in_place_wahbm(char *in,
                                 uint32_t query_value,
                                 uint32_t *R,
                                 uint32_t num_records,
                                 int time,
                                 int quiet,
                                 char *bim)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R = gt_records_compressed_in_place_wahbm(wf,
                                                                 R,
                                                                 num_records,
                                                                 query_value,
                                                                 &wf_R);
    uint32_t *ints;
    uint32_t len_ints = compressed_in_place_wah_to_ints(wf_R,
                                                            len_wf_R,
                                                            &ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields, bim);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

int gt_wahbm(char *in,
             uint32_t query_value,
             uint32_t *R,
             uint32_t num_records,
             int time,
             int quiet,
             char *bim)

{
    start();
    struct wah_file wf = init_wahbm_file(in);
    uint32_t *wf_R;
    uint32_t len_wf_R = gt_records_wahbm(wf,
                                             R,
                                             num_records,
                                             query_value,
                                             &wf_R);

    uint32_t *ints;
    uint32_t len_ints = wah_to_ints(wf_R,len_wf_R,&ints);
    stop();
    if (time != 0)
        fprintf(stderr,"%lu\n", report());

    if (quiet == 0)
        print_result(len_ints, ints, wf.num_fields, bim);

    free(ints);
    free(wf_R);
    fclose(wf.file);

    return 0;
}

void print_result(uint32_t len,
                  uint32_t *R,
                  uint32_t num_fields,
                  char *bim)
{

	/* OLD WAY *********************
    uint32_t i,j, bit_i = 0;
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


    NEW WAY */

    uint32_t i,j,line_idx,bytes, bit_i = 0;

	struct quick_file_info qfile;
	struct output_buffer outbuf;


	init_out_buf(&outbuf, NULL);
	quick_file_init(bim, &qfile);

	for (i=0; i < len; ++i) {
		bytes = R[i];
		if (bytes == 0) continue; /* skip a bunch of ops if you can */
		for (j=0; j < 32; j++) {
			if (bytes & 1 << (31 - j)) {
				line_idx = i*32+j;
				append_out_buf(&outbuf, qfile.lines[line_idx], qfile.line_lens[line_idx]);
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
