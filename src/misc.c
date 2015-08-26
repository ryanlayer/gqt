#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "bcf.h"
#include "ubin.h"
#include "wahbm.h"

int misc_help();
int plt_ubin(char *in, char *out);
int plt_vcf(char *in, char *out);
int ubin_plt(char *in, char *out);
int ubin_wahbm(char *in, char *out, char *full_cmd);
int ubin_wahbm16(char *in, char *out);
int ubin_wah(char *in, char *out);
int vcf_plt(char *in,
            char *out,
            uint32_t num_fields,
            uint32_t num_records);
int plt_invert(char *in, char *out);
int plt_invert_ubin(char *in, char *out);
int wahbm_pca(char *in, char *out);
int wahbm_hamm(char *in, char *out);
int wahbm_shared(char *in, char *out);
//int top_n_matches(char *in, char *db, uint32_t num_matches);
int top_n_matches(char *in, uint32_t num_matches);
int speed_check(char *in);

int misc(int argc, char **argv, char *full_cmd)
{
    if (argc < 2) return misc_help();

    int c;
    char *in = NULL, *out = NULL, *bim = NULL, *vid = NULL, *db = NULL;
    uint32_t num_fields = 0, num_records = 0, num_matches = 0;
    int i_is_set = 0, 
        o_is_set = 0, 
        d_is_set = 0, 
        f_is_set = 0, 
        b_is_set = 0, 
        v_is_set = 0, 
        n_is_set = 0, 
        r_is_set = 0; 

    while((c = getopt (argc, argv, "hi:o:f:r:b:v:n:d:")) != -1) {
        switch (c) {
            case 'd':
                d_is_set = 1;
                db = optarg;
                break;
            case 'n':
                n_is_set = 1;
                num_matches = atoi(optarg);
                break;
            case 'v':
                v_is_set = 1;
                vid = optarg;
                break;
            case 'b':
                b_is_set = 1;
                bim = optarg;
                break;
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
                misc_help();
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
                misc_help();
                return 1;
        }
    }

    char *type = argv[0];



    if (strcmp(type, "top-n") == 0) {
        //if ( (i_is_set == 0) || (n_is_set ==0) || (d_is_set == 0) {
        if ( (i_is_set == 0) || (n_is_set ==0) ) {
            //printf("top-n requires an input file, n\n");
            printf("top-n requires an input file and n to be set\n");
            return misc_help();
        }

        return top_n_matches(in, num_matches);
    }


    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return misc_help();
    } 
    if (o_is_set == 0) {
        printf("Output file is not set\n");
        return misc_help();
    }
    if (strcmp(type, "vcf-plt") == 0) {
        if (f_is_set == 0) {
            printf("Number of fields is not set\n");
            return misc_help();
        }
        if (r_is_set == 0) {
            printf("Number of records is not set\n");
            return misc_help();
        }

        return vcf_plt(in, out, num_fields, num_records);
    } 
    if (strcmp(type, "rotate-hlubin") == 0) {
        if (f_is_set == 0) {
            printf("Number of fields is not set\n");
            return misc_help();
        }
        if (r_is_set == 0) {
            printf("Number of records is not set\n");
            return misc_help();
        }

        rotate_gt(num_fields,
                  num_records,
                  in,
                  out);
        return 0;
    }

    if (strcmp(type, "plt-ubin") == 0)  return plt_ubin(in, out);
    if (strcmp(type, "plt-vcf") == 0)  return plt_vcf(in, out);
    if (strcmp(type, "plt-invert") == 0)  return plt_invert(in, out);
    if (strcmp(type, "plt-invert-ubin") == 0)  return plt_invert_ubin(in, out);
    if (strcmp(type, "ubin-plt") == 0) return ubin_plt(in, out);
    if (strcmp(type, "ubin-wahbm") == 0) return ubin_wahbm(in, out, full_cmd);
    if (strcmp(type, "ubin-wahbm16") == 0) return ubin_wahbm16(in, out);
    if (strcmp(type, "ubin-wah") == 0) return ubin_wah(in, out);
    if (strcmp(type, "pca") == 0) return wahbm_pca(in, out);
    if (strcmp(type, "hamm") == 0) return wahbm_hamm(in, out);
    if (strcmp(type, "shared") == 0) return wahbm_shared(in, out);
    if (strcmp(type, "top-n") == 0) return wahbm_pca(in, out);
    if (strcmp(type, "speed-check") == 0) return speed_check(in);

    return 1;
}

int misc_help()
{
    printf("usage:   gqt covert <type> -i <input file> -o <output file>\n"
           "         plt-invert        Switch records to fields\n"
           "         plt-invert-ubin   Switch records to fields\n"
           "         plt-ubin          Plain text to uncompress binary\n"
           "         plt-vcf           Plain text to VCF\n"
           "         ubin-plt          Uncompressed binary to plain text\n"
           "         ubin-wahbm        Uncompressed binary to WAH bitmap\n"
           "         ubin-wahbm16      Uncompressed binary to 16-bit WAH "
                                       "bitmap\n"
           "         ubin-wah          Uncompressed binary to WAH \n"
           "         vcf-plt           VCF to by-variant plain text\n"
           "         bcf-wahbm         BCF to by-individual sorted WAH bitmap\n"
           "         ped-db            PED to SQLite3 database\n"
           "         rotate-hlubin     Rotate a ubin w/o a header to one w/\n"
           "         pca               Run PCA/\n"
           "         hamm              Get pair-wise hamming distance/\n"
           "         shared            Get pair-wise sharing size/\n"
           "         speed-check       Test ops/\n"
           "         top-n             Find the top n matching hets/\n"
           "         -v                VID output file name"
                                       "(required for bcf-wahbm)\n"
           "         -b                BIM output file name"
                                       "(required for bcf-wahbm)\n"
           "         -r                Number of records "
                                       "(required for vcf-plt and bcf-wahbm)\n"
           "         -f                Number of fields "
                                       "(required for vcf-plt and bcf-wahbm)\n"
    );

    return 0;
}

int ubin_wahbm16(char *in, char *out)
{
    return convert_file_by_name_ubin_to_wahbm16(in, out);
}

int ubin_wahbm(char *in, char *out, char *full_cmd)
{
    return convert_file_by_name_ubin_to_wahbm(in, out, full_cmd);
}

int ubin_plt(char *in, char *out)
{
    return convert_file_by_name_ubin_to_plt(in, out);
}

int ubin_wah(char *in, char *out)
{
    return convert_file_by_name_ubin_to_wah(in, out);
}

int plt_ubin(char *in, char *out)
{
    return convert_file_by_name_plt_to_ubin(in, out);
}

int plt_vcf(char *in, char *out)
{
    return convert_file_by_name_plt_to_vcf(in, out);
}

int plt_invert(char *in, char *out)
{
    return convert_file_by_name_invert_plt(in, out);
}

int plt_invert_ubin(char *in, char *out)
{
    return convert_file_by_name_invert_plt_to_ubin(in, out);
}


int vcf_plt(char *in,
            char *out,
            uint32_t num_fields,
            uint32_t num_records)
{
    return convert_file_by_name_vcf_to_plt(in, num_fields, num_records, out);
}

int wahbm_pca(char *in,
              char *out)
{
    return wahbm_pca_by_name(in, out);
}

int wahbm_hamm(char *in,
               char *out)
{
    return wahbm_hamm_dist_by_name(in, out);
}

int wahbm_shared(char *in,
               char *out)
{
    return wahbm_shared_by_name(in, out);
}



int speed_check(char *in)
{
    return wahbm_speed_check(in);
}

int top_n_matches(char *in, uint32_t num_matches)
{
    return wahbm_top_n_matches_by_name(in, num_matches);
}


