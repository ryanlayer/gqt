#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage()
{
    fprintf(stderr,
        "usage:\tor_gt_records_from_uint <options>\n"
        "\t\t-u\tInput variant or individual uncompressed interger file name\n"
        "\t\t-R\tNumber of records\n"
        "\t\t-F\tNumber of fields\n"
        "\t\t-n\tNumber of records to OR\n"
        "\t\t-r\tCSV record IDs\n"
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name;
    char *record_ids;
    int num_r, num_records, num_fields;
    int u_is_set = 0,
        R_is_set = 0, 
        F_is_set = 0, 
        r_is_set = 0, 
        n_is_set = 0; 

    while ((c = getopt (argc, argv, "u:r:n:R:F:")) != -1) {
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
                num_r = atoi(optarg);
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

    if (  (r_is_set == 0) || 
          (u_is_set == 0) || 
          (R_is_set == 0) || 
          (F_is_set == 0) || 
          (n_is_set == 0)) {
        usage();
        return 1;
    }

    int R[num_r];
    parse_cmd_line_int_csv(R, num_r, record_ids);

    struct uint_file u_file = init_uint_file(in_file_name,
                                             num_records,
                                             num_fields);

    unsigned int *r;
    or_uint_records(u_file, num_r, R, &r);

    fclose(u_file.file);
    free(u_file.line);

    /*
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);

    int i,j,k=0;
    for (i = 0; i < num_ints_per_record; ++i) {
        unsigned int c = r[i];
        for (j = 0; j < 16; ++j) {
            if (k == u_file.num_fields)
                break;

            int gt = (c >> (30 - (j*2))) & 3;

            printf("%d ", gt);

            k+=1;
        }
    }
    printf("\n");
    */
}
