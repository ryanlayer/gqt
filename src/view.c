#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

int view_help();
int view_plt(char *in);
int view_ubin(char *in, int record_id);
int view_wahbm(char *in);

int view(int argc, char **argv)
{
    if (argc < 2) return view_help();

    int c;
    char *in, *out;
    int i_is_set = 0,
        r_is_set; 
    int record_id = -1;

    while ((c = getopt (argc, argv, "hi:r:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in = optarg;
                break;
            case 'r':
                r_is_set = 1;
                record_id = atoi(optarg);
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
                return view_help();
            default:
                return view_help();
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        fprintf(stderr, "Input file is not set\n");
        return view_help();
    } 

    if (strcmp(type, "plt") == 0)  return view_plt(in);
    else if (strcmp(type, "ubin") == 0) return view_ubin(in, record_id);
    else if (strcmp(type, "wahbm") == 0) return view_wahbm(in);

    return 1;
}

int view_help()
{
    fprintf(stderr,
            "%s %s\n"
            "usage:   gqt view <type> -i <input file>\n"
            "         plt   Plain text\n"
            "         ubin  Uncompressed binary\n"
            "         wahbm WAH-encoded bitmap\n"
            "         -r    Record number to print\n",
            PROGRAM_NAME, VERSION);

    return 1;
}

int view_plt(char *in)
{
    uint32_t num_printed = print_by_name_plt(in, NULL, 0);
    return 0;
}
int view_ubin(char *in, int record_id)
{
    if (record_id != -1) {
        uint32_t R[1] = { record_id };
        uint32_t num_printed = print_by_name_ubin(in, R, 1, 0);
    } else {
        uint32_t num_printed = print_by_name_ubin(in, NULL, 0, 0);
    }
    return 0;
}

int view_wahbm(char *in)
{
    uint32_t num_printed = print_by_name_wahbm(in, NULL, 0, 0);
    return 0;
}
