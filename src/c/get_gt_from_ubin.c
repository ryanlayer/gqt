#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage()
{
    fprintf(stderr,
        "usage:\tint_to_ubin <options>\n"
        "\t\t-u\tInput variant or individual uncompressed binary file name\n"
        "\t\t-n\tNumber of genotypes (record/filed pairs) to get\n"
        "\t\t-r\tCSV record IDs\n"
        "\t\t-f\tCSV field IDs\n"
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name;
    char *record_ids, *field_ids;
    int num_gt;
    int u_is_set = 0,
        r_is_set = 0, 
        n_is_set = 0, 
        f_is_set = 0; 

    while ((c = getopt (argc, argv, "u:r:f:n:")) != -1) {
        switch (c) {
            case 'n':
                n_is_set = 1;
                num_gt = atoi(optarg);
                break;
            case 'f':
                f_is_set = 1;
                field_ids = optarg;
                break;
            case 'r':
                r_is_set = 1;
                record_ids = optarg;
                break;
            case 'u':
                u_is_set = 1;
                in_file_name = optarg;
                break;
            case 'h':
                usage();
                return 1;
            case '?':
                if ( (optopt == 'f') || (optopt == 'r') ||
                     (optopt == 'u') )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                usage();
                return 1;
        }
    }

    if (  (f_is_set == 0) || 
          (r_is_set == 0) || 
          (u_is_set == 0) || 
          (n_is_set == 0)) {
        usage();
        return 1;
    }

    int R[num_gt];
    parse_cmd_line_int_csv(R, num_gt, record_ids);

    int F[num_gt];
    parse_cmd_line_int_csv(F, num_gt, field_ids);
    
    /*
    FILE *i_file = fopen(in_file_name, "rb");
    if (!i_file) {
        fprintf(stderr, "Unable to open %s\n",in_file_name);
        return 1;
    }
    */
    struct ubin_file u_file = init_ubin_file(in_file_name);

    int gts[num_gt];

    get_ubin_genotypes(u_file, num_gt, R, F, gts);

    int i;
    for (i = 0; i < num_gt; ++i)
        printf("%d %d %d\n", R[i], F[i], gts[i]); 

    fclose(u_file.file);
}
