#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <sysexits.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

#include "bcf.h"
#include "ped.h"
#include "genotq.h"

int convert_help();
int bcf_wahbm_metadata(char *bcf_file_name,
                 char *wah_out,
                 char *bim_out,
                 char *vid_out,
                 char *tmp_dir,
                 uint32_t num_fields,
                 uint32_t num_records,
                 char *full_cmd);
int bcf_wahbm_offset(char *bcf_file_name,
                     char *wah_out,
                     char *offset_out,
                     char *vid_out,
                     char *tmp_dir,
                     uint32_t num_fields,
                     uint32_t num_records,
                     char *full_cmd);

int ped_ped(char *bcf_file_name,
            char *ped_file_name,
            uint32_t col,
            char *ped_db_file_name);

int convert(int argc, char **argv, char *full_cmd)
{
    if (argc < 2) return convert_help();

    int c;
    char *bcf_file_name=NULL,
         *gqt_file_name=NULL,
         *ped_db_file_name=NULL,
         *bim_file_name=NULL,
         *off_file_name=NULL,
         *vid_file_name=NULL,
         *tmp_dir=NULL,
         *ped_file_name=NULL;

    uint32_t num_fields = 0, num_records = 0, col = 2;
    int i_is_set = 0, 
        p_is_set = 0, 
        c_is_set = 0, 
        t_is_set = 0, 
        r_is_set = 0, 
        f_is_set = 0, 
        G_is_set = 0, 
        V_is_set = 0, 
        O_is_set = 0, 
        B_is_set = 0, 
        D_is_set = 0;

    while((c = getopt (argc, argv, "i:p:c:t:f:r:G:V:O:B:D:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                bcf_file_name = optarg;
                break;
            case 'p':
                p_is_set = 1;
                ped_file_name = optarg;
                break;
            case 'c':
                col = atoi(optarg);
                break;
            case 't':
                t_is_set = 1;
                tmp_dir = optarg;
                break;
            case 'f':
                f_is_set = 1;
                num_fields = atoi(optarg);
                break;
            case 'r':
                r_is_set = 1;
                num_records = atoi(optarg);
                break;
            case 'G':
                G_is_set = 1;
                gqt_file_name = optarg;
                break;
            case 'V':
                V_is_set = 1;
                vid_file_name = optarg;
                break;
            case 'O':
                O_is_set = 1;
                off_file_name = optarg;
                break;
            case 'B':
                B_is_set = 1;
                bim_file_name = optarg;
                break;
            case 'D':
                D_is_set = 1;
                ped_db_file_name = optarg;
                break;
            case 'h':
                convert_help();
                return 1;
            case '?':
                if ( (optopt == 'i') || 
                     (optopt == 'p') ||
                     (optopt == 'c') ||
                     (optopt == 't') ||
                     (optopt == 'f') ||
                     (optopt == 'r') ||
                     (optopt == 'G') ||
                     (optopt == 'V') ||
                     (optopt == 'O') ||
                     (optopt == 'B') ||
                     (optopt == 'D') )
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

    // Make sure that the BCF file is specified and accessible 
    if (i_is_set == 0) {
        fprintf(stderr, "BCF/VCF/VCF.GZ file is not set\n");
        return convert_help();
    } else {
        if ( access( bcf_file_name, F_OK) == -1 )
            err(EX_NOINPUT,
                "Error accessing BCF/VCF/VCF.GZ file '%s'",
                bcf_file_name);
    }

    if (strcmp(type, "bcf") == 0) {
        // If the user does not specify fields and rows then look for index
        if ( (f_is_set == 0) || (r_is_set == 0) ) {

            fprintf(stderr,"Attempting to autodetect number of variants "
                    "and samples from %s\n", bcf_file_name);
            //Try and auto detect the sizes, need the index
            tbx_t *tbx = NULL;
            hts_idx_t *idx = NULL;
            htsFile *fp    = hts_open(bcf_file_name,"rb");
            if ( !fp ) {
                err(EX_DATAERR, "Could not read file: %s", bcf_file_name);
            }

            bcf_hdr_t *hdr = bcf_hdr_read(fp);
            if ( !hdr ) {
                err(EX_DATAERR, "Could not read the header: %s", bcf_file_name);
            }

            if (hts_get_format(fp)->format==vcf) {
                tbx = tbx_index_load(bcf_file_name);
                if ( !tbx ) { 
                    err(EX_NOINPUT,
                        "Could not load index: %s.csi",
                        bcf_file_name);
                }
            } else if ( hts_get_format(fp)->format==bcf ) {
                idx = bcf_index_load(bcf_file_name);
                if ( !idx ) {
                    err(EX_NOINPUT,
                        "Could not load index: %s.csi",
                        bcf_file_name);
                }
            } else {
                err(EX_NOINPUT,
                    "Could not detect the file type as VCF or BCF: %s",
                    bcf_file_name);
            }

            num_fields = hdr->n[BCF_DT_SAMPLE];

            num_records = 0;
            const char **seq;
            int nseq;
            seq = tbx ? tbx_seqnames(tbx, &nseq) : 
                    bcf_index_seqnames(idx, hdr, &nseq);
            int i;
            uint32_t sum = 0;
            for (i = 0; i < nseq; ++i) {
                uint64_t records, v;
                hts_idx_get_stat(tbx ? tbx->idx: idx, i, &records, &v);
                num_records += records;
            }

            fprintf(stderr, "Number of variants:%u\tNumber of samples:%u\n",
                    num_records, num_fields);
            free(seq);
            hts_close(fp);
            bcf_hdr_destroy(hdr);
            if (idx)
                hts_idx_destroy(idx);
            if (tbx)
                tbx_destroy(tbx);
        }

        if (G_is_set == 0) {
            gqt_file_name  = (char*)malloc(strlen(bcf_file_name) + 5);
            if (!gqt_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(gqt_file_name,bcf_file_name);
            strcat(gqt_file_name, ".gqt");
        }

        if (B_is_set == 0) {
            bim_file_name  = (char*)malloc(strlen(bcf_file_name) + 5);
            if (!bim_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(bim_file_name,bcf_file_name);
            strcat(bim_file_name, ".bim");
        }

        if (O_is_set == 0) {
            off_file_name  = (char*)malloc(strlen(bcf_file_name) + 5); 
            if (!off_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(off_file_name, bcf_file_name);
            strcat(off_file_name, ".off");
        }

        if (V_is_set == 0) {
            vid_file_name  = (char*)malloc(strlen(bcf_file_name) + 5);
            if (!vid_file_name)
                err(EX_OSERR, "malloc error");
            strcpy(vid_file_name, bcf_file_name);
            strcat(vid_file_name, ".vid");
        }

        if (t_is_set == 0) {
            tmp_dir  = (char*)malloc(3*sizeof(char)); // "./\0"
            if (!tmp_dir  )
                err(EX_OSERR, "malloc error");
            strcpy(tmp_dir,"./");
        }

        return convert_file_by_name_bcf_to_wahbm_metadata_offset(
                bcf_file_name,
                num_fields,
                num_records,
                gqt_file_name,
                bim_file_name,
                off_file_name,
                vid_file_name,
                tmp_dir,
                full_cmd);
    } 

    if (strcmp(type, "ped") == 0)  {
        if (D_is_set == 0) {
            if (p_is_set == 1) {
                ped_db_file_name  = (char*)malloc(strlen(ped_file_name) + 4);
                if (!ped_db_file_name)
                    err(EX_OSERR, "malloc error");
                strcpy(ped_db_file_name, ped_file_name);
                strcat(ped_db_file_name, ".db");
            } else {
                ped_db_file_name = (char*)malloc(strlen(bcf_file_name) + 4);
                if (!ped_db_file_name)
                    err(EX_OSERR, "malloc error");
                strcpy(ped_db_file_name,bcf_file_name);
                strcat(ped_db_file_name, ".db");
            }
      }

      return ped_ped(bcf_file_name, ped_file_name, col, ped_db_file_name);
    }
    return convert_help();
}

int convert_help()
{
    fprintf(stderr,
            "%s v%s\n"
            "usage:   gqt convert <type> -i <input VCF/VCF.GZ/BCF file>\n"
            "     types:\n"
            "         bcf         create a GQT index\n"
            "         ped         create sample phenotype database\n\n"
            "     options:\n"
            "         -p           PED file name (opt.)\n"
            "         -c           Sample name column in PED (Default 2)\n"
            "         -G           GQT output file name (opt.)\n"
            "         -V           VID output file name (opt.)\n"
            "         -O           OFF output file name (opt.)\n"
            "         -B           BIM output file name (opt.)\n"
            "         -D           PED DB output file name (opt.)\n"
            "         -r           Number of variants (opt. with index)\n"
            "         -f           Number of samples (opt. with index)\n"
            "         -t           Tmp working directory(./ by default)\n",
            PROGRAM_NAME, VERSION);

    return EX_USAGE;
}


int ped_ped(char *bcf_file_name,
            char *ped_file_name,
            uint32_t col,
            char *ped_db_file_name)
{
    return convert_file_by_name_ped_to_db(bcf_file_name,
                                          ped_file_name,
                                          col,
                                          ped_db_file_name);
}

int bcf_wahbm_metadata(char *bcf_file_name,
                       char *wah_out,
                       char *bim_out,
                       char *vid_out,
                       char *tmp_dir,
                       uint32_t num_fields,
                       uint32_t num_records,
                       char *full_cmd)
{
    return convert_file_by_name_bcf_to_wahbm_bim(bcf_file_name,
                                                 num_fields,
                                                 num_records,
                                                 wah_out,
                                                 bim_out,
                                                 vid_out,
                                                 tmp_dir,
                                                 full_cmd);
}

int bcf_wahbm_offset(char *bcf_file_name,
                     char *wah_out,
                     char *offset_out,
                     char *vid_out,
                     char *tmp_dir,
                     uint32_t num_fields,
                     uint32_t num_records,
                     char *full_cmd)
{
    return convert_file_by_name_bcf_to_wahbm_offset(bcf_file_name,
                                                    num_fields,
                                                    num_records,
                                                    wah_out,
                                                    offset_out,
                                                    vid_out,
                                                    tmp_dir,
                                                    full_cmd);
}
