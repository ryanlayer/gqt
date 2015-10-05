#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sysexits.h>
#include <stdlib.h>

#include "off.h"
#include "genotq.h"

struct off_file *new_off_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples)
{
    struct off_file *o = (struct off_file *) malloc(sizeof(struct off_file));
    if (!o)
        err(EX_OSERR, "malloc error");

    o->gqt_header = new_gqt_file_header('o', 
                                        full_cmd,
                                        num_variants,
                                        num_samples);

    o->file_name = strdup(file_name);

    o->file = fopen(file_name,"wb");
    if (!o->file)
        err(EX_CANTCREAT, "Cannot create OFF file '%s'", file_name);

    if (fwrite(o->gqt_header,
               sizeof(struct gqt_file_header),
               1,
               o->file) != 1)
        err(EX_IOERR, "Error writing header to OFF file '%s'", file_name);

    return o;
}

void destroy_off_file(struct off_file *o)
{
    if (fclose(o->file))
        err(EX_IOERR, "Error closeing OFF file '%s'", o->file_name); 

    free(o->gqt_header);
    free(o->file_name);
    o = NULL;
}

struct off_file *open_off_file(char *file_name)
{
    struct off_file *o = (struct off_file *) malloc(sizeof(struct off_file));
    if (!o)
        err(EX_OSERR, "malloc error");

    o->file_name = strdup(file_name);
    o->file = fopen(file_name,"rb+");

    if (!(o->file))
        err(EX_NOINPUT, "Cannot open OFF file '%s'", file_name);

    o->gqt_header = read_gqt_file_header(o->file_name, o->file);

    if ( !((o->gqt_header->marker[0] == 'G') &&
           (o->gqt_header->marker[1] == 'Q') && 
           (o->gqt_header->marker[2] == 'T')) )
        errx(EX_NOINPUT, "File '%s' is not a GQT file.", file_name);

    if (o->gqt_header->type != 'o')
        errx(EX_NOINPUT, "File '%s' is not a OFF file.", file_name);

    o->offsets = (uint64_t *)
        malloc((o->gqt_header->num_variants)*sizeof(uint64_t));
    if (!(o->offsets))
        err(EX_OSERR, "malloc error");


    size_t fr = fread(o->offsets,
                      sizeof(uint64_t),
                      o->gqt_header->num_variants,
                      o->file);
    check_file_read(o->file_name, o->file, o->gqt_header->num_variants, fr);

    return o;
}

void add_to_off_file(struct off_file *o, uint64_t next_off)
{
    if (fwrite(&next_off, sizeof(uint64_t), 1, o->file) != 1)
        err(EX_IOERR, "Error writing header to OFF file '%s'", o->file_name);
}
