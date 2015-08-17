#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sysexits.h>

#include "vid.h"

//{{{struct vid_file *new_vid_file(char *file_name,
struct vid_file *new_vid_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples)
{
    struct vid_file *v = (struct vid_file *) malloc(sizeof(struct vid_file));

    v->file_name = strdup(file_name);

    v->file = fopen(file_name,"wb");
    if (!v->file)
        err(EX_CANTCREAT, "Cannot create VID file \"%s\"", file_name);

    v->gqt_header = new_gqt_file_header('v', 
                                        full_cmd,
                                        num_variants,
                                        num_samples);

    v->vids = NULL;

    if (fwrite(v->gqt_header,
               sizeof(struct gqt_file_header),
               1,
               v->file) != 1)
        err(EX_IOERR, "Error writing header to VID file \"%s\"", file_name);

    return v;
}
//}}}

//{{{void write_vid(struct vid_file *v, uint32_t vid)
void write_vid(struct vid_file *v, uint32_t vid)
{
    if (fwrite(&vid, sizeof(uint32_t), 1, v->file) != 1)
        err(EX_IOERR, "Error writing to \"%s\"", v->file_name); 
}
//}}}

//{{{struct vid_file *open_vid_file(char *file_name)
struct vid_file *open_vid_file(char *file_name)
{
    struct vid_file *v = (struct vid_file *) malloc(sizeof(struct vid_file));

    v->file_name = strdup(file_name);
    v->file = fopen(file_name,"rb");

    if (!(v->file))
        err(EX_NOINPUT, "Cannot open VID file \"%s\"", file_name);

    v->gqt_header = read_gqt_file_header(v->file_name, v->file);

    if ( !((v->gqt_header->marker[0] == 'G') &&
           (v->gqt_header->marker[1] == 'Q') && 
           (v->gqt_header->marker[2] == 'T')) )
        errx(EX_NOINPUT, "File '%s' is not a GQT file.", file_name);

    if (v->gqt_header->type != 'v')
        errx(EX_NOINPUT, "File '%s' is not a VID file.", file_name);


    v->vids = NULL;

    return v;
}
//}}}

//{{{void load_vid_data(struct vid_file *v)
void load_vid_data(struct vid_file *v)
{
    if (v->vids != NULL)
        errx(EX_SOFTWARE, 
             "VID data has already been loaded for file '%s'.", v->file_name);

    v->vids = (uint32_t *) malloc(v->gqt_header->num_variants*sizeof(uint32_t));
    if (!v->vids)
        err(EX_OSERR, "malloc error");

    if (fseek(v->file, sizeof(struct gqt_file_header), SEEK_SET))
        err(EX_IOERR, "Error seeking to data in VID file '%s'.", v->file_name);

    size_t fr = fread(v->vids,
                      sizeof(uint32_t),
                      v->gqt_header->num_variants,
                      v->file);
    check_file_read(v->file_name, v->file, v->gqt_header->num_variants, fr);
}
//}}}

//{{{void destroy_vid_file(struct vid_file *v)
void destroy_vid_file(struct vid_file *v)
{
    if (fclose(v->file))
        err(EX_IOERR, "Error closeing VID file '%s'", v->file_name); 

    free(v->file_name);

    if (v->vids != NULL) 
        free(v->vids);

    free(v);
    v = NULL;
}
//}}}
