#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage()
{
    fprintf(stderr,
        "usage:\tor_gt_fields_from_uint <options>\n"
        "\t\t-u\tInput variant or individual uncompressed interger file name\n"
        "\t\t-R\tNumber of records\n"
        "\t\t-F\tNumber of fields\n"
        "\t\t-n\tNumber of fields to OR\n"
        "\t\t-f\tCSV field IDs\n"
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name;
    char *field_ids;
    int num_f, num_records, num_fields;
    int u_is_set = 0,
        R_is_set = 0, 
        F_is_set = 0, 
        f_is_set = 0, 
        n_is_set = 0; 

    while ((c = getopt (argc, argv, "u:f:n:R:F:")) != -1) {
        switch (c) {
            case 'F':
                F_is_set = 1;
                num_fields = atoi(optarg);
                break;
            case 'R':
                R_is_set = 1;
                num_records = atoi(optarg);
                break;
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
                if ( (optopt == 'f') || (optopt == 'f') ||
                     (optopt == 'F') || (optopt == 'R') ||
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
          (R_is_set == 0) || 
          (F_is_set == 0) || 
          (n_is_set == 0)) {
        usage();
        return 1;
    }

    int F[num_f];
    parse_cmd_line_int_csv(F, num_f, field_ids);

    struct uint_file u_file = init_uint_file(in_file_name,
                                             num_records,
                                             num_fields);

    unsigned int *f;
    or_uint_fields(u_file, num_f, F, &f);

    fclose(u_file.file);

    int i;
    for (i = 0; i < u_file.num_records; ++i)
            printf("%u ", f[i]);
    printf("\n");
}
