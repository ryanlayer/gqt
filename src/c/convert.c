#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

int convert_help();
int plt_ubin(char *in, char *out);
int ubin_wahbm(char *in, char *out);
int ubin_wahbm16(char *in, char *out);
int ubin_wah(char *in, char *out);
int vcf_plt(char *in,
            char *out,
            unsigned int num_fields,
            unsigned int num_records);

int convert(int argc, char **argv)
{
    if (argc < 2) return convert_help();

    int c;
    char *in, *out;
    unsigned int num_fields, num_records;
    int i_is_set = 0, 
        o_is_set = 0, 
        f_is_set = 0, 
        r_is_set = 0; 

    while ((c = getopt (argc, argv, "hi:o:f:r:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in = optarg;
                break;
            case 'o':
                o_is_set = 1;
                out = optarg;
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
                convert_help();
                return 1;
            case '?':
                if ( (optopt == 'i') || 
                     (optopt == 'f') ||
                     (optopt == 'r') ||
                     (optopt == 'o') )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                convert_help();
                return 1;
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return convert_help();
    } 
    if (o_is_set == 0) {
        printf("Output file is not set\n");
        return convert_help();
    }
    if (strcmp(type, "vcf-plt") == 0) {
        if (f_is_set == 0) {
            printf("Number of fields is not set\n");
            return convert_help();
        }
        if (r_is_set == 0) {
            printf("Number of records is not set\n");
            return convert_help();
        }

        return vcf_plt(in, out, num_fields, num_records);
    } 
    if (strcmp(type, "plt-ubin") == 0)  return plt_ubin(in, out);
    if (strcmp(type, "ubin-wahbm") == 0) return ubin_wahbm(in, out);
    if (strcmp(type, "ubin-wahbm16") == 0) return ubin_wahbm16(in, out);
    if (strcmp(type, "ubin-wah") == 0) return ubin_wah(in, out);

    return 1;
}

int convert_help()
{
    printf("usage:   gtq covert <type> -i <input file> -o <output file>\n"
           "         plt-ubin     Plain text to uncompress binary\n"
           "         ubin-wahbm   Uncompressed binary to WAH bitmap\n"
           "         ubin-wahbm16 Uncompressed binary to 16-bit WAH bitmap\n"
           "         ubin-wah     Uncompressed binary to WAH \n"
           "         vcf-plt      VCF to by-variant plain text\n"
           "         -r           Number of records (required for vcf-plt)\n"
           "         -f           Number of fields (required for vcf-plt)\n"
    );

    return 0;
}

int ubin_wahbm16(char *in, char *out)
{
    return convert_file_by_name_ubin_to_wahbm16(in, out);
}

int ubin_wahbm(char *in, char *out)
{
    return convert_file_by_name_ubin_to_wahbm(in, out);
}

int ubin_wah(char *in, char *out)
{
    return convert_file_by_name_ubin_to_wah(in, out);
}

int plt_ubin(char *in, char *out)
{
    return convert_file_by_name_plt_to_ubin(in, out);
}

int vcf_plt(char *in,
            char *out,
            unsigned int num_fields,
            unsigned int num_records)
{
    return convert_file_by_name_vcf_to_plt(in, num_fields, num_records, out);
}
