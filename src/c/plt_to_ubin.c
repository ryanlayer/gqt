#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\%s <options>\n"
        "\t\t-i\tInput variant or individual plain text file name\n"
        "\t\t-o\tOuput file name\n", prog
    );
}

int main(int argc, char **argv)
{
    char *progr = argv[0];
    int c;
    char *in_file_name, *out_file_name;
    int i_is_set = 0, 
        o_is_set = 0; 

    while ((c = getopt (argc, argv, "hi:o:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in_file_name = optarg;
                break;
            case 'o':
                o_is_set = 1;
                out_file_name = optarg;
                break;
            case 'h':
                usage(prog);
                return 1;
            case '?':
                if ( (optopt == 'i') || (optopt == 'o') )
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

    if (  (o_is_set == 0) || (i_is_set == 0) ) {
        usage(prog);
        return 1;
    }

    plt_to_ubin(in_file_name, out_file_name);
    struct plt_file plt_f = init_plt_file(in_file_name);
    
    
#if 0    
    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First value is the number of fields per record
    fwrite(&(plt_f.num_fields), sizeof(int), 1, o_file);

    // Second value is the number of records
    fwrite(&(plot_f.num_records), sizeof(int), 1, o_file);

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
#endif
}

