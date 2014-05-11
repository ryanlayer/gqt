#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genotq.h"

static int test_1 = 0;
static char *line;
static int line_len;
static FILE *file;
static int read_state = 0;

//{{{ void init_int_genotype_reader(char *file_name, int num_gt)
void init_int_genotype_reader(char *file_name, int num_gt)
{

    line_len = num_gt*2 + 1;
    line = (char *) malloc(line_len*sizeof(char));
    file = fopen(file_name, "r");
    read_state = 0;
}
//}}}

//{{{void destroy_int_genotype_reader()
void destroy_int_genotype_reader()
{
    fclose(file);
    free(line);
}
//}}}

//{{{ int get_next_int_genotype(int *line_num, int *gt_num, int *gt)
int get_next_int_genotype(int *line_num, int *gt_num, int *gt)
{
    // read state: 0 get a line
    // read state: 1 get first char
    // read state: 2 get next char
    // read state: 3 EOF

    char *pch;

    if (read_state == 0) {
        char *r = fgets(line, line_len, file);
        if (r == NULL) {
            read_state = 3;
            return 1;
        } else {
            read_state = 1;
        }
    }

    if (read_state == 1) { 
        pch = strtok(line, " ");
       *gt_num = 0;

        read_state = 2;
    } else if (read_state == 2) {
        pch = strtok(NULL, " ");
       *gt_num = *gt_num + 1;
    }

    if (pch != NULL) {
        *gt = atoi(pch);
        return 0;
    } else {
        read_state = 0;
        *line_num = *line_num + 1;
        return get_next_int_genotype(line_num, gt_num, gt);
    }
}
//}}}

struct uint_file init_uint_file(char *file_name, 
                                int num_records,
                                int num_fields)
{
    struct uint_file uf;
    uf.line_len = num_fields*2 + 1;
    uf.line = (char *) malloc(line_len*sizeof(char));
    uf.file = fopen(file_name, "rb");

    return uf;
}

int or_uint_records(struct uint_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r)
{
    return 0;
}


//{{{ struct ubin_file init_ubin_file(char *file_name)
struct ubin_file init_ubin_file(char *file_name)
{
    struct ubin_file uf;

    uf.file = fopen(file_name, "rb");

    if (!uf.file) {
        fprintf(stderr, "Unable to open %s\n", file_name);
        return uf;
    }

    // Jump to the begining of the file to grab the record size
    fseek(uf.file, 0, SEEK_SET);
    fread(&uf.num_fields,sizeof(unsigned int),1,uf.file);
    fread(&uf.num_records,sizeof(unsigned int),1,uf.file);

    return uf;

}
///}}}

//{{{ int get_ubin_genotypes(struct ubin_file u_file,
int get_ubin_genotypes(struct ubin_file u_file,
                       int num_gts,
                       int *record_ids,
                       int *field_ids,
                       int *gts)
{
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    int i;
    int field_int_id, field_offset;
    unsigned int target_int;
    for (i = 0; i < num_gts; ++i) {
        field_int_id = field_ids[i]/16;
        field_offset = field_ids[i] - field_int_id*16;

        unsigned int target_int;
        // skip to the int that contains the target field
        fseek(u_file.file, 8 + // skip the record and field size field
                    record_ids[i]*num_bytes_per_record + // skip to the reccord
                    field_int_id*4,// skip to the propper int 
                    SEEK_SET);
        fread(&target_int,sizeof(unsigned int),1,u_file.file);

        // shift the int over to get the target field
        gts[i] = (target_int >> (30 - field_offset*2)) & 3;
    }

    return 0;
}
//}}}

//{{{ int or_ubin_fields(struct ubin_file u_file, 
int or_ubin_fields(struct ubin_file u_file, 
                   int num_f,
                   int *field_ids,
                   unsigned int **r)
{


    int i,j;
    *r = (unsigned int *) malloc(u_file.num_records*sizeof(unsigned int));
    int *record_ids = (int *) malloc(num_f*sizeof(int));
    int *gts = (int *) malloc(num_f*sizeof(int));
    for (i = 0; i < u_file.num_records; ++i) {
        for (j = 0; j < num_f; ++j) {
            record_ids[j] = i;
        }

        get_ubin_genotypes(u_file,
                           num_f,
                           record_ids,
                           field_ids,
                           gts);


        (*r)[i]=0;
        for (j = 0; j < num_f; ++j) {
            (*r)[i]= (*r)[i] | gts[j];
        }
    }

    free(record_ids);

    return 0;
}
//}}}

//{{{ int or_ubin_records(struct ubin_file u_file, 
int or_ubin_records(struct ubin_file u_file, 
                    int num_r,
                    int *record_ids,
                    unsigned int **r)
{
    int num_ints_per_record = 1 + ((u_file.num_fields - 1) / 16);
    int num_bytes_per_record = num_ints_per_record * 4;

    // create enough space in r, and set all the values to zero
    *r = (unsigned int *) calloc(num_ints_per_record, sizeof(unsigned int));

    unsigned int *c = (unsigned int *)
        malloc(num_ints_per_record*sizeof(unsigned int));

    int i,j;
    for (i = 0; i < num_r; ++i) {

        // skip to the target record and read in the full record
        fseek(u_file.file, 8 + // skip the record and field size field
                    record_ids[i]*num_bytes_per_record, // skip to the reccord
                    SEEK_SET);
        fread(c,sizeof(unsigned int),num_bytes_per_record,u_file.file);

        //and it
        for (j = 0; j < num_ints_per_record; ++j) {
            (*r)[j] = (*r)[j] | c[j];
        }
    }

    return 0;
}
//}}}

//{{{ void parse_cmd_line_int_csv(int *I,
void parse_cmd_line_int_csv(int *I,
                            int num_I,
                            char *cmd_line_arg)
{
    char *pch;
    pch = strtok(cmd_line_arg,",");
    int i;
    for (i = 0; i < num_I; ++i){
        I[i] = atoi(pch);
        pch = strtok(NULL,",");
    }
}
//}}}
