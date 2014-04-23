#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage()
{
    fprintf(stderr,
        "usage:\tor_gt_fields_from_ubin <options>\n"
        "\t\t-u\tInput variant or individual uncompressed binary file name\n"
        "\t\t-n\tNumber of fields to OR\n"
        "\t\t-f\tCSV field IDs\n"
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name;
    char *field_ids;
    int num_f;
    int u_is_set = 0,
        f_is_set = 0, 
        n_is_set = 0; 

    while ((c = getopt (argc, argv, "u:f:n:")) != -1) {
        switch (c) {
            case 'n':
                n_is_set = 1;
                num_f = atoi(optarg);
                break;
            case 'f':
                f_is_set = 1;
                field_ids = optarg;
                break;
            case 'u':
                u_is_set = 1;
                in_file_name = optarg;
                break;
            case 'h':
                usage();
                return 1;
            case '?':
                if ( (optopt == 'f') || (optopt == 'n') ||
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
          (u_is_set == 0) || 
          (n_is_set == 0)) {
        usage();
        return 1;
    }

    int F[num_f];
    parse_cmd_line_int_csv(F, num_f, field_ids);

    struct ubin_file u_file = init_ubin_file(in_file_name);

    unsigned int *r;
    or_ubin_fields(u_file, num_f, F, &r);

    fclose(u_file.file);

    int i;
    for (i = 0; i < u_file.num_records; ++i) {
        if (i !=0)
            printf(" ");
        printf("%u",r[i]);

    }
    printf("\n");
}
