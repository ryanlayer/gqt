#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage()
{
    fprintf(stderr,
        "usage:\tint_to_ubin <options>\n"
        "\t\t-i\tInput variant or individual integer-encoded file name\n"
        "\t\t-r\tNumber of records (rows)\n"
        "\t\t-f\tNumber of fields (columns) per record\n"
        "\t\t-o\tOuput file name\n"
    );
}

int main(int argc, char **argv)
{
    int c;
    char *in_file_name, *out_file_name;
    int num_fields, num_records;
    int f_is_set = 0,
        i_is_set = 0, 
        r_is_set = 0, 
        o_is_set = 0; 

    while ((c = getopt (argc, argv, "f:i:r:o:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in_file_name = optarg;
                break;
            case 'o':
                o_is_set = 1;
                out_file_name = optarg;
                break;
            case 'f':
                f_is_set = 1;
                num_fields = atoi(optarg);
                break;
            case 'r':
                r_is_set = 1;
                num_records = atoi(optarg);
                break;
            case 'h':
                usage();
                return 1;
            case '?':
                if ( (optopt == 'f') || (optopt == 'o') ||
                     (optopt == 'i') || (optopt == 'r') )
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

    if (  (f_is_set == 0) || (i_is_set == 0) ||
        (r_is_set == 0) || (o_is_set == 0)) {
        usage();
        return 1;
    }
        
    
    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First value is the number of fields per record
    fwrite(&num_fields, sizeof(int), 1, o_file);

    // Second value is the number of records
    fwrite(&num_records, sizeof(int), 1, o_file);

    init_int_genotype_reader(in_file_name, num_fields);
    int line_num = 0, gt_num = 0, gt = 0;

    int i = 1;

    unsigned int curr = 0;
    int last_line_num = 0;
    while (get_next_int_genotype(&line_num, &gt_num, &gt) == 0) {

        if (line_num != last_line_num) {
            if ( i > 1) 
                //printf("%u\n", curr);
                fwrite(&curr, sizeof(unsigned int), 1, o_file);
            //printf("\n");
            curr = 0;
            i = 1;
        }

        curr += gt << (32 - i*2);

        ++i;

        if (i == 16) {
            fwrite(&curr, sizeof(unsigned int), 1, o_file);
            //printf("%u\n", curr);
            curr = 0;
            i = 1;
        }

        last_line_num = line_num;
    }

    if (i > 1)
        fwrite(&curr, sizeof(unsigned int), 1, o_file);

    destroy_int_genotype_reader();
    fclose(o_file);
}

