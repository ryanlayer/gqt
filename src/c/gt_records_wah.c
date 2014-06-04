#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\t%s<options>\n"
        "\t\t-i\tInput variant or individual uncompressed binary file name\n"
        "\t\t-n\tNumber of records to OR\n"
        "\t\t-r\tCSV record IDs\n",
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
        n_is_set = 0; 

    while ((c = getopt (argc, argv, "i:r:n:")) != -1) {
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

    struct wah_file wf = init_wah_file(in_file_name);

    unsigned int *wf_R;
    unsigned int len_wf_R = gt_records_wah(wf, R, num_r, 0, &wf_R);

    return 0;
}
