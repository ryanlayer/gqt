#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

int view_help();
int view_plt(char *in);
int view_ubin(char *in);
int view_wah(char *in);
int view_wahbm(char *in);

int view(int argc, char **argv)
{
    if (argc < 2) return view_help();

    int c;
    char *in, *out;
    int i_is_set = 0; 

    while ((c = getopt (argc, argv, "hi:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in = optarg;
                break;
            case 'h':
                view_help();
                return 1;
            case '?':
                if (optopt == 'i') 
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                view_help();
                return 1;
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return view_help();
    } 

    if (strcmp(type, "plt") == 0)  return view_plt(in);
    else if (strcmp(type, "ubin") == 0) return view_ubin(in);
    else if (strcmp(type, "wah") == 0) return view_wah(in);
    else if (strcmp(type, "wahbm") == 0) return view_wahbm(in);

    return 1;
}

int view_help()
{
    printf("usage:   gtq view <type> -i <input file>\n"
           "         plt-ubin    Plain text to uncompress binary\n"
           "         ubin-wahbm  Uncompressed binary to WAH bitmap\n"
           "         ubin-wah    Uncompressed binary to WAH \n"
    );

    return 0;
}

int view_plt(char *in)
{
    unsigned int num_printed = print_by_name_plt(in, NULL, 0);
    return 0;
}
int view_ubin(char *in)
{
    unsigned int num_printed = print_by_name_ubin(in, NULL, 0, 0);
    return 0;
}
int view_wah(char *in)
{
    unsigned int num_printed = print_by_name_wah(in, NULL, 0, 0);
    return 0;
}
int view_wahbm(char *in)
{
    unsigned int num_printed = print_by_name_wahbm(in, NULL, 0, 0);
    return 0;
}
