/**
 * @file genotq.c
 * @Author Ryan Layer (ryan.layer@gmail.com)
 * @date May, 2014
 * @brief Functions for converting and opperation on various encoding
 * strategies
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <math.h>
#include <limits.h>
#include "genotq.h"
#include "plt.h"



//{{{ struct plt_file init_plt_file(char *file_name)
struct plt_file init_plt_file(char *file_name)
{
    struct plt_file plt_f;

    plt_f.file = fopen(file_name, "r");

    if (!plt_f.file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        abort();
    }

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    int gt;

    read = getline(&line, &len, plt_f.file);

    plt_f.num_fields = atoi(line);

    read = getline(&line, &len, plt_f.file);

    plt_f.num_records = atoi(line);

    plt_f.header_offset = ftell(plt_f.file);

    free(line);

    return plt_f;
}
///}}}

//{{{int convert_file_by_name_vcf_to_plt(char *in_file_name, char
int convert_file_by_name_vcf_to_plt(char *in_file_name,
                                    unsigned int num_fields,
                                    unsigned int num_records,
                                    char *out_file_name)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *pch;
    int gt;

    FILE *i_file = fopen(in_file_name, "r");
    if (!i_file) {
        fprintf(stderr, "Unable to open %s\n",in_file_name);
        return 1;
    }

    FILE *o_file = fopen(out_file_name, "w");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }


    fprintf(o_file, "%u\n%u\n", num_fields, num_records);

    while ((read = getline(&line, &len, i_file)) != -1) {
        if (line[0] != '#') {

            // Skip the first 9 fields
            
            int i;
            pch = strtok(line, "\t");
            for  (i =0; i < 8; ++i) {
                pch = strtok(NULL, "\t");
            }
            
            // get the first genotype
            pch = strtok(NULL, "\t");
            i = 0;
            while (pch != NULL) {
                if ((pch[0] == '0') && (pch[2] == '0'))
                    gt = 0;
                else if ((pch[0] == '1') && (pch[2] == '0'))
                    gt = 1;
                else if ((pch[0] == '0') && (pch[2] == '1'))
                    gt = 1;
                else if ((pch[0] == '1') && (pch[2] == '1'))
                    gt = 2;
                else
                    gt = 3;

                if (i != 0)
                    fprintf(o_file, " ");

                fprintf(o_file, "%d", gt);
                pch = strtok(NULL, "\t");
                i+=1;
            }
            fprintf(o_file, "\n");
        }
    }

    free(line);
    fclose(i_file);
    fclose(o_file);
    return 0;
}
//}}}

//{{{ void convert_file_by_name_plt_to_vcf(char *in_file_name, 
int convert_file_by_name_plt_to_vcf(char *in_file_name, char *out_file_name)
{
    struct plt_file pf = init_plt_file(in_file_name);
    int r = convert_file_plt_to_vcf(pf, out_file_name);
    fclose(pf.file);
    return r;
}

int convert_file_plt_to_vcf(struct plt_file pf, char *out_file_name)
{

    char *line = NULL;
    size_t len;
    ssize_t read;
    unsigned int *packed_ints;

    FILE *o_file = fopen(out_file_name, "w");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First print the VCF headers
    fprintf(o_file,
            "##fileformat=VCFv4.1\n"
            "##INFO=<ID=N,Type=String>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##contig=<ID=1,length=3000000000>\n");

    // Print field headers
    fprintf(o_file,
            "#CHROM\t"
            "POS\t"
            "ID\t"
            "REF\t"
            "ALT\t"
            "QUAL\t"
            "FILTER\t"
            "INFO\t"
            "FORMAT\t");

    
    unsigned int i,j,g;
    for (i = 0; i < pf.num_fields; ++i){
        if (i != 0)
            fprintf(o_file,"\t");
        fprintf(o_file,"I%u", i);
    }
    fprintf(o_file,"\n");
    
    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);
    for (i = 0; i < pf.num_records; ++i) {
        fprintf(o_file,
            "1\t"
            "%u\t"
            "V%u\t"
            "A\t"
            "T\t"
            "100\t"
            "PASS\t"
            "N=A\t"
            "GT\t", i+1, i+1);

        read = getline(&line, &len, pf.file);
        for (j = 0; j < pf.num_fields; ++j) {
            if (j != 0)
                fprintf(o_file,"\t");
            g = ((int)line[j*2] - 48);
            if (g == 0)
                fprintf(o_file,"0|0");
            else if (g == 1)
                fprintf(o_file,"0|1");
            else if (g == 2)
                fprintf(o_file,"1|1");

        }
        fprintf(o_file,"\n");
    }

    fclose(o_file);

    return 0;
}
//}}}

//{{{ void convert_file_by_name_plt_to_ubin(char *in_file_name, 
int convert_file_by_name_plt_to_ubin(char *in_file_name, char *out_file_name)
{
    struct plt_file pf = init_plt_file(in_file_name);
    int r = convert_file_plt_to_ubin(pf, out_file_name);
    fclose(pf.file);
    return r;
}

int convert_file_plt_to_ubin(struct plt_file pf, char *out_file_name)
{

    char *line = NULL;
    unsigned int i,j,count;
    size_t len;
    ssize_t read;
    int g[16];
    unsigned int *packed_ints;

    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    // First value is the number of fields per record
    fwrite(&(pf.num_fields), sizeof(int), 1, o_file);

    // Second value is the number of records
    fwrite(&(pf.num_records), sizeof(int), 1, o_file);


    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);
    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        //fprintf(stderr, "read:%lu\n", read);
        //fprintf(stderr, "%s", line);
        unsigned int num_packed_ints  = plt_line_to_packed_ints(line,
                                                                pf.num_fields,
                                                                &packed_ints);

        for (j = 0; j < num_packed_ints; ++j)
            fwrite(&(packed_ints[j]), sizeof(unsigned int), 1, o_file);

        free(packed_ints);
    }

    fclose(o_file);

    return 0;
}
//}}}

//{{{ unsigned int plt_line_to_packed_ints(char *line,
unsigned int plt_line_to_packed_ints(char *line,
                                     unsigned int num_fields, 
                                     unsigned int **packed_ints)
{
    unsigned int len = 1 + ((num_fields - 1) / 16);
    *packed_ints = (unsigned *) malloc((len)*sizeof(unsigned int));
    int i, two_bit_count, pack_int_count = 0;

    //printf("%s\n", line);
    
    int g[16];

    memset(g,0,sizeof(g));
    two_bit_count = 0;
    for (i = 0; i < num_fields; ++i) {
        g[i%16] = ((int)line[i*2] - 48);
        //printf("%d %d %c\n", i, i*2, line[i*2]);

        two_bit_count +=1;

        if (two_bit_count == 16) {
            int j;
            /*
            for (j = 0; j < two_bit_count; ++j) 
                printf("%d ", g[j]);
            printf("\n");
            */
            (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);
            pack_int_count += 1;
            memset(g,0,sizeof(g));
            two_bit_count = 0;
        }
        
    }


    /*
    for (i = 0; i < two_bit_count; ++i) 
        printf("%d ", g[i]);
    printf("\n");
    */

    if (two_bit_count > 0)
        (*packed_ints)[pack_int_count] = pack_2_bit_ints(g, 16);

    return len;
}
//}}}

//{{{ int convert_file_by_name_invert_plt(char *in_file_name, char
int convert_file_by_name_invert_plt(char *in_file_name, char *out_file_name)
{
    struct plt_file pf = init_plt_file(in_file_name);
    FILE *o_pf = fopen(out_file_name, "w");

    // get the path to the ouput to write temp files
    char full_path[PATH_MAX+1];
    char *ptr = realpath(in_file_name, full_path);
    
    char orig_path[strlen(ptr)];
    char *o = strcpy(orig_path, full_path);

    char *pch, *last;
    pch = strtok(ptr, "/");
    while (pch != NULL) {
        last = pch;
        pch = strtok(NULL, "/");
    }

    // the number of temp file will be equal to the number of records (the old
    // number of fields), to make sure there is enough room to hold the string
    // name
    unsigned int max_digit_len = (unsigned int)log10((double)pf.num_fields)+1;

    char temp_path[strlen(orig_path) - strlen(last) + 4 + max_digit_len + 1];
    unsigned int path_len = strlen(orig_path) - strlen(last);
    pch = strncpy(temp_path, orig_path, path_len);
    pch = strncpy(temp_path + path_len, ".tmp", 4);
    path_len += 4;


    fprintf(o_pf, "%u\n", pf.num_records);
    fprintf(o_pf, "%u\n", pf.num_fields);

     

    char *line = NULL;
    ssize_t read;
    size_t len;
    unsigned int i,j,g,n;

    // clear out any old files
    for (j = 0; j < pf.num_fields; ++j) {
        n = sprintf(temp_path + path_len, "%u", j);
        FILE *tmp_file = fopen(temp_path, "w");
        fclose(tmp_file);
    }

    // put each line into a seperate file
    fseek(pf.file, pf.header_offset, SEEK_SET);
    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        for (j = 0; j < pf.num_fields; ++j) {
            n = sprintf(temp_path + path_len, "%u", j);
            FILE *tmp_file = fopen(temp_path, "a");
            if (i!= 0)
                fprintf(tmp_file," ");
            g = ((int)line[j*2] - 48);
            fprintf(tmp_file,"%u", g);
            fclose(tmp_file);
        }
    }

    // combine and remove tmp files
    for (j = 0; j < pf.num_fields; ++j) {
        n = sprintf(temp_path + path_len, "%u", j);
        FILE *tmp_file = fopen(temp_path, "r");
        read = getline(&line, &len, tmp_file);
        fprintf(o_pf,"%s\n", line);
        fclose(tmp_file);
        remove(temp_path);
    } 

        
    fclose(pf.file);
    fclose(o_pf);
    return 0;
}
//}}}

//{{{ int convert_file_by_name_invert_plt_to_ubin(char *in_file_name,
int convert_file_by_name_invert_plt_to_ubin(char *in_file_name,
                                            char *out_file_name)
{
    struct plt_file pf = init_plt_file(in_file_name);

    FILE *o_file = fopen(out_file_name, "wb");
    if (!o_file) {
        fprintf(stderr, "Unable to open %s\n",out_file_name);
        return 1;
    }

    unsigned int n_num_fields, n_num_records;
    unsigned int i, two_bit_i = 0;
    unsigned int *ubin = NULL;
    char *line = NULL;
    size_t len;
    ssize_t read;
    
    // put each line into a seperate file
    fseek(pf.file, pf.header_offset, SEEK_SET);
    for (i = 0; i < pf.num_records; ++i) {
        fprintf(stderr, "I:%u\n", i);
        read = getline(&line, &len, pf.file);
        two_bit_i = invert_plt_to_ubin(line,
                                       pf.num_fields,
                                       pf.num_records,
                                       &n_num_fields,
                                       &n_num_records,
                                       two_bit_i,
                                       &ubin);
    }


    // First value is the number of fields per record
    fwrite(&(n_num_fields), sizeof(unsigned int), 1, o_file);

    // Second value is the number of records
    fwrite(&(n_num_records), sizeof(unsigned int), 1, o_file);

    unsigned int num_ints_per_record = 1 + (((n_num_fields) - 1) / 16);
    for (i = 0; i < n_num_records; ++i) {
        fprintf(stderr, "O:%u\n", i);
        fwrite(&(ubin[i*num_ints_per_record]),
               sizeof(unsigned int),
               num_ints_per_record,
               o_file);
    }
        
    fclose(pf.file);
    fclose(o_file);

    return 0;
}
//}}}

//{{{ unsigned int invert_plt_to_ubin(char *line,
unsigned int invert_plt_to_ubin(char *line,
                                unsigned int num_fields,
                                unsigned int num_records,
                                unsigned int *new_num_fields,
                                unsigned int *new_num_records,
                                unsigned int field_i,
                                unsigned int **ubin)
{

    *new_num_fields = num_records;
    *new_num_records = num_fields;
    unsigned int num_ints_per_record = 1 + (((*new_num_fields) - 1) / 16);
    unsigned int int_i = field_i / 16;
    unsigned int two_bit_i = field_i % 16;

    if (*ubin == NULL)
        *ubin = (unsigned int *) calloc(num_ints_per_record*(*new_num_records), 
                                        sizeof(unsigned int));
    
    unsigned int i,g;
    for (i = 0; i < num_fields; ++i) {
        g = ((int)line[i*2] - 48);
        /*
        fprintf(stderr, "i:%u\t"
                        "int_i:%u\t"
                        "two_bit_i:%u\t"
                        "offset:%u\t"
                        "\n",
                        i, int_i, two_bit_i,
                        i*num_ints_per_record+int_i
                        );
        */
        (*ubin)[i*num_ints_per_record + int_i] += g << (30 - two_bit_i*2);
        //fprintf(stderr, "u:%u\n", (*ubin)[i*num_ints_per_record + int_i]);
    }
    
    return field_i + 1;
}
//}}}

//{{{unsigned int pack_2_bit_ints(int *ints, int num_ints)
unsigned int pack_2_bit_ints(int *ints, int num_ints)
{
    int i;
    unsigned int r = 0;
    for (i = 0; i < num_ints; ++i) {
       r = r | ((ints[i] & 3) << (30 - i*2));
    }
    
    return r;
}
//}}}

//{{{ unsigned int get_plt_record(struct plt_file pf,
unsigned int get_plt_record(struct plt_file pf,
                            unsigned int plt_record,
                            unsigned int **plt)
{

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *plt = (unsigned int *) calloc(num_ints_per_record, sizeof(unsigned int));

    long line_len = pf.num_fields*2*sizeof(char);

    fseek(pf.file, pf.header_offset + line_len*plt_record, SEEK_SET);
    read = getline(&line, &len, pf.file);

    unsigned int plt_size = plt_line_to_packed_ints(line,
                                                    pf.num_fields, 
                                                    plt);
    free(line);
    return plt_size;
}
//}}}

//{{{ unsigned int count_range_records_plt(struct plt_file pf,
unsigned int count_range_records_plt(struct plt_file pf,
                                     unsigned int *record_ids,
                                     unsigned int num_r,
                                     unsigned int start_test_value,
                                     unsigned int end_test_value,
                                     unsigned int **R)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,val;

    *R = (unsigned int *) calloc(pf.num_fields,sizeof(unsigned int));

    long line_len = pf.num_fields*2*sizeof(char);
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        for (j = 0; j < pf.num_fields; ++j) {
            val = (unsigned int)line[j*2] - 48;            
            if ( (val > start_test_value) && (val < end_test_value) )
                (*R)[j] += 1;
        }
    }

    free(line);

    return pf.num_fields;
}
//}}}

//{{{ unsigned int range_records_plt(struct plt_file pf,
unsigned int range_records_plt(struct plt_file pf,
                              unsigned int *record_ids,
                              unsigned int num_r,
                              unsigned int start_test_value,
                              unsigned int end_test_value,
                              unsigned int **R)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *R = (unsigned int *) malloc(num_ints_per_record*sizeof(unsigned int));

    for (i = 0; i < num_ints_per_record; ++i)
        (*R)[i] = -1;

    unsigned int int_i, bit_i;


    long line_len = pf.num_fields*2*sizeof(char);
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        int_i = 0;
        bit_i = 0;
        for (j = 0; j < pf.num_fields; ++j) {
            // clear the bit if fails criteria
            unsigned int val = (unsigned int)line[j*2] - 48;            
            if (!(val >= start_test_value && val < end_test_value))
                (*R)[int_i] = (*R)[int_i] & ~(1 << (31 - bit_i));

            bit_i += 1;
            if (bit_i == 32) {
                int_i += 1;
                bit_i = 0;
            }
        }
    }

    free(line);

    return num_ints_per_record;
}
//}}}

//{{{ START eq, gt, ne, gte, lt, lte: records_plt
//{{{ unsigned int eq_records_plt(struct plt_file pf,
unsigned int eq_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value, test_value+1, R);
}
//}}}

//{{{ unsigned int gt_count_records_plt(struct plt_file pf,
unsigned int gt_count_records_plt(struct plt_file pf,
                                  unsigned int *record_ids,
                                  unsigned int num_r,
                                  unsigned int test_value,
                                  unsigned int **R)
{
    // TODO: need constants for upper bound.
    return count_range_records_plt(pf,
                                   record_ids,
                                   num_r,
                                   test_value,
                                   4,
                                   R);
}
//}}}

//{{{ unsigned int ne_records_plt(struct plt_file pf,
unsigned int ne_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;

    unsigned int num_ints_per_record = 1 + ((pf.num_fields - 1) / 32);

    *R = (unsigned int *) malloc(num_ints_per_record*sizeof(unsigned int));

    for (i = 0; i < num_ints_per_record; ++i)
        (*R)[i] = -1;

    unsigned int int_i, bit_i;


    long line_len = pf.num_fields*2*sizeof(char);
    for (i = 0; i < num_r; ++i) {
        fseek(pf.file, pf.header_offset + line_len*record_ids[i], SEEK_SET);

        read = getline(&line, &len, pf.file);

        int_i = 0;
        bit_i = 0;
        for (j = 0; j < pf.num_fields; ++j) {
            // clear the bit
            unsigned int val = ((unsigned int)line[j*2] - 48);
            if  (!(val != test_value))
                (*R)[int_i] = (*R)[int_i] & ~(1 << (31 - bit_i));

            bit_i += 1;
            if (bit_i == 32) {
                int_i += 1;
                bit_i = 0;
            }
        }
    }

    free(line);

    return num_ints_per_record;
}
//}}}

//{{{ unsigned int gt_records_plt(struct plt_file pf,
unsigned int gt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value+1, 4, R);
}
//}}}

//{{{ unsigned int gte_records_plt(struct plt_file pf,
unsigned int gte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, test_value, 4, R);
}
//}}}

//{{{ unsigned int lt_records_plt(struct plt_file pf,
unsigned int lt_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, 0, test_value, R);
}
//}}}

//{{{ unsigned int lte_records_plt(struct plt_file pf,
unsigned int lte_records_plt(struct plt_file pf,
                            unsigned int *record_ids,
                            unsigned int num_r,
                            unsigned int test_value,
                            unsigned int **R)
{
    // TODO: need constants for upper bound.
    return range_records_plt(pf, record_ids, num_r, 0, test_value+1, R);
}
//}}}
//}}} END eq, gt, ne, gte, lt, lte: records_plt

//{{{ unsigned int print_plt(struct plt_file pf,
unsigned int print_plt(struct plt_file pf,
                       unsigned int *record_ids,
                       unsigned int num_r)
{
    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    unsigned int i,j,bit;
    unsigned int num_printed = 0;


    long line_len = pf.num_fields*2*sizeof(char);

    if ( (num_r == 0) || record_ids == NULL ) {
        while ((read = getline(&line, &len, pf.file)) != -1) {
            printf("%s", line);
            num_printed += 1;
        }
    } else {
        for (i = 0; i < num_r; ++i) {
            fseek(pf.file,
                  pf.header_offset + line_len*record_ids[i],
                  SEEK_SET);

            read = getline(&line, &len, pf.file);
            printf("%s", line);
            num_printed += 1;
        }
    }

    free(line);

    return num_printed;
}
//}}}

//{{{unsigned int print_by_name_plt(char *pf_file_name,
unsigned int print_by_name_plt(char *pf_file_name,
                               unsigned int *record_ids,
                               unsigned int num_r)
{
    struct plt_file pf = init_plt_file(pf_file_name);
    return print_plt(pf, record_ids, num_r); 

}
//}}}

//{{{ unsigned int plt_to_bitmap_wah(unsigned int *U,
unsigned int plt_to_bitmap_wah(char *plt,
                               unsigned int plt_len,
                               unsigned int **W,
                               unsigned int **wah_sizes)
{

    unsigned int *ubin;
    unsigned int ubin_len = plt_line_to_packed_ints(plt, plt_len, &ubin);

    unsigned int wah_len = ubin_to_bitmap_wah(ubin,
                                              ubin_len,
                                              plt_len, // one bit per field
                                              W,
                                              wah_sizes);
    free(ubin);

    return wah_len;
}
//}}}

