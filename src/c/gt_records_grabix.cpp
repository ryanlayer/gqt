#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"
#include "grabix.h"

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\t%s<options>\n"
        "\t\t-i\tInput variant or individual bgzip and grabix indexed file\n"
        "\t\t-n\tNumber of records to OR\n"
        "\t\t-r\tCSV record IDs\n"
        "\t\t-q\tQuiet mode\n",
        prog
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name;
    char *record_ids;
    char *prog = argv[0];
    int num_r;
    int i_is_set = 0,
        r_is_set = 0, 
        n_is_set = 0, 
        q_is_set = 0; 

    while ((c = getopt (argc, argv, "i:r:n:q")) != -1) {
        switch (c) {
            case 'n':
                n_is_set = 1;
                num_r = atoi(optarg);
                break;
            case 'r':
                r_is_set = 1;
                record_ids = optarg;
                break;
            case 'i':
                i_is_set = 1;
                in_file_name = optarg;
                break;
            case 'q':
                q_is_set = 1;
                break;
            case 'h':
                usage(prog);
                return 1;
            case '?':
                if ( (optopt == 'f') || (optopt == 'r') ||
                     (optopt == 'i') )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                usage(prog);
                return 1;
        }
    }

    if (  (r_is_set == 0) || 
          (i_is_set == 0) || 
          (n_is_set == 0)) {
        usage(prog);
        return 1;
    }

    unsigned int R[num_r];
    parse_cmd_line_int_csv(R, num_r, record_ids);

    index_info index;
    BGZF *bgzf_fp;

    if (bgzf_init(in_file_name, index, &bgzf_fp)) {
        cerr << "[grabix] error loading index/opending file" << endl;
        return 1;
    }

    char *line;
    size_t len;
    len = get_line(index, bgzf_fp, 1, &line);

    unsigned int num_fields = atoi(line);

    free(line);

    len = get_line(index, bgzf_fp, 2, &line);

    unsigned int num_records = atoi(line);

    free(line);



    unsigned int num_ints_per_record = 1 + ((num_fields - 1) / 32);

    unsigned int *G =
            (unsigned int *) malloc(num_ints_per_record*sizeof(unsigned int));

    int i, j;
    for (i = 0; i < num_ints_per_record; ++i)
        G[i] = -1;

    unsigned int int_i, bit_i;
    for (i = 0; i < num_r; ++i) {
        len = get_line(index, bgzf_fp, R[i]+3, &line);
        
        //cout << line << endl;

        int_i = 0;
        bit_i = 0;
        for (j = 0; j < num_fields; ++j) {
            //cout << ((unsigned int)line[j*2] - 48) << endl;
            // clear the bit
            if  ( !( 0 < ((unsigned int)line[j*2] - 48))) 
                G[int_i] = G[int_i] & ~(1 << (31 - bit_i));

            bit_i += 1;
            if (bit_i == 32) {
                int_i += 1;
                bit_i = 0;
            }
        }
        free(line);
    }

    bgzf_close(bgzf_fp);

    if (q_is_set == 0) {
        for(i = 0; i < num_ints_per_record; ++i) {
            if (i != 0)
                printf(" ");
            printf("%u", G[i]);
        }
        printf("\n");
    }


    return 0;
}
