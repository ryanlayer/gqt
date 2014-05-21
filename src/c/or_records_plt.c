#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\t%s <options>\n"
        "\t\t-i\tInput variant or individual plain text file name\n"
        "\t\t-n\tNumber of records to OR\n"
        "\t\t-r\tCSV record IDs\n", 
        prog
    );
}

int main(int argc, char **argv)
{
    int c;
    char *prog = argv[0];
    char *in_file_name;
    char *record_ids;
    int num_r, num_records, num_fields;
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

    int R[num_r];
    parse_cmd_line_int_csv(R, num_r, record_ids);

    struct plt_file pf = init_plt_file(in_file_name);

    int *G = (int *) calloc(pf.num_fields, sizeof(int));

    int r = or_records_plt(pf, R, num_r, G);

    int i;
    for (i = 0; i < pf.num_fields; ++i) {
        if (i != 0)
            printf(" ");
        printf("%d", G[i]);
    }
    printf("\n");

    fclose(pf.file);
    free(G);

    return 0;
}
