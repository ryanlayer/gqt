#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

int convert_help();
int bcf_wahbm(char *in,
              char *wah_out,
              char *bim_out,
              char *vid_out,
              char *tmp_dir,
              uint32_t num_fields,
              uint32_t num_records);
int ped_db(char *in, char *out);

int convert(int argc, char **argv)
{
    if (argc < 2) return convert_help();

    int c;
    char *in=NULL, *out=NULL, *bim=NULL, *vid=NULL, *tmp_dir=NULL;
    uint32_t num_fields, num_records;
    int i_is_set = 0, 
        o_is_set = 0, 
        f_is_set = 0, 
        b_is_set = 0, 
        v_is_set = 0, 
        t_is_set = 0, 
        r_is_set = 0; 

    while((c = getopt (argc, argv, "hi:o:f:r:b:v:t:")) != -1) {
        switch (c) {
            case 't':
                t_is_set = 1;
                tmp_dir = optarg;
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
                convert_help();
                return 1;
            case '?':
                if ( (optopt == 'i') || 
                     (optopt == 'f') ||
                     (optopt == 'r') ||
                     (optopt == 't') ||
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

    if (strcmp(type, "bcf") == 0) {
        if (f_is_set == 0) {
            // TODO auto detect from BCF index
            printf("Number of fields is not set\n");
            return convert_help();
        }
        if (r_is_set == 0) {
            // TODO auto detect from BCF index
            printf("Number of records is not set\n");
            return convert_help();
        }
        if (o_is_set == 0) {
            out  = (char*)malloc(strlen(in) + 5); // 5 for ext and \0
            strcpy(out,in);
            strcat(out, ".gqt");
        }
        if (b_is_set == 0) {
            bim  = (char*)malloc(strlen(in) + 5); // 5 for ext and \0
            strcpy(bim,in);
            strcat(bim, ".bim");
        }
        if (v_is_set == 0) {
            vid  = (char*)malloc(strlen(in) + 5); // 5 for ext and \0
            strcpy(vid,in);
            strcat(vid, ".vid");
        }
        if (t_is_set == 0) {
            tmp_dir  = (char*)malloc(3*sizeof(char)); // "./\0"
            strcpy(tmp_dir,"./");
        }

        int r = bcf_wahbm(in, out, bim, vid, tmp_dir, num_fields, num_records);

        /*
        if (vid != NULL)
            free(vid);
        if (bim != NULL)
            free(bim);
        if (out != NULL)
            free(out);
        */


        return r;
    } 

    if (strcmp(type, "ped") == 0)  {
      if (o_is_set == 0) {
            out  = (char*)malloc(strlen(in) + 4); // 4 for ext and \0
            strcpy(out,in);
            strcat(out, ".db");
      }
      return ped_db(in, out);
    
    }
    return convert_help();

}

int convert_help()
{
    printf("usage:   gqt convert <type> -i <input file>\n"
           "     types:\n"
           "         bcf         create a WAH index of a VCF/VCF.gz/BCF\n"
           "         ped         create SQLite3 database of a PED file\n\n"
           "     options:\n"
           "         -o           Output file name\n"
           "         -v           VID output file name (opt.)\n"
           "         -b           BIM output file name (opt.)\n"
           "         -r           Number of variants (req. for bcf)\n"
           "         -f           Number of samples (req. for bcf)\n"
           "         -t           Tmp working directory(./ by defualt)\n"
           );

    return 0;
}


int ped_db(char *in, char *out)
{
    return convert_file_by_name_ped_to_db(in, out);
}

int bcf_wahbm(char *in,
              char *wah_out,
              char *bim_out,
              char *vid_out,
              char *tmp_dir,
              uint32_t num_fields,
              uint32_t num_records)
{
    return convert_file_by_name_bcf_to_wahbm_bim(in,
                                                 num_fields,
                                                 num_records,
                                                 wah_out,
                                                 bim_out,
                                                 vid_out,
                                                 tmp_dir);
}
