#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\%s <options>\n"
        "\t\t-i\tInput variant or individual uncompressed binary file\n"
        "\t\t-o\tOuput file name\n", prog
    );
}

int main(int argc, char **argv)
{
    char *prog = argv[0];
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

    int r = convert_file_by_name_ubin_to_wahbm(in_file_name, out_file_name);
    
    if (r == 0)
        return 0;
    else {
        fprintf(stderr, "Error creating file\n"); 
        return 1;
    }
}
