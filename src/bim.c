#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sysexits.h>
#include <assert.h>
#include <htslib/knetfile.h>

#include "bim.h"

 /* 
 * The file is :
 * uncompressed size     ( sizeof(uint64_t))
 * compressed size       ( sizeof(uint64_t))
 * header size           ( sizeof(uint64_t))
 * md line lengths       ( bcf_f->num_records*sizeof(uint64_t))
 * compressed data 
 */
struct bim_file_header *new_bim_file_header(uint64_t u_size,
                                            uint64_t c_size,
                                            uint64_t h_size,
                                            uint64_t *md_line_lens)
{
    struct bim_file_header *h = 
            (struct bim_file_header *) malloc(sizeof(struct bim_file_header));

    h->u_size = u_size;
    h->c_size = c_size;
    h->h_size = h_size;
    h->md_line_lens = md_line_lens;

    return h;
}

struct bim_file *new_bim_file(char *file_name,
                              char *full_cmd,
                              uint32_t num_variants,
                              uint32_t num_samples,
                              uint64_t u_size,
                              uint64_t c_size,
                              uint64_t h_size,
                              uint64_t *md_line_lens)
{
    struct bim_file *b = (struct bim_file *) malloc(sizeof(struct bim_file));

    b->gqt_header = new_gqt_file_header('b', 
                                        full_cmd,
                                        num_variants,
                                        num_samples);

    b->bim_header = new_bim_file_header(u_size,
                                        c_size,
                                        h_size,
                                        md_line_lens);

    b->file_name = strdup(file_name);

    b->type = BIM_LOCAL;

    b->file.local = fopen(file_name,"wb");
    if (!b->file.local)
        err(EX_CANTCREAT, "Cannot create BIM file \"%s\"", file_name);

    if (fwrite(b->gqt_header,
               sizeof(struct gqt_file_header),
               1,
               b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&u_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&c_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&h_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(md_line_lens,
               sizeof(uint64_t),
               num_variants,
               b->file.local) != num_variants)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    b->data_start = ftell(b->file.local);

    return b;
}

void destroy_bim_file(struct bim_file *b)
{
    if (b->type == BIM_LOCAL) {
        if (fclose(b->file.local))
            err(EX_IOERR, "Error closeing BIM file '%s'", b->file_name); 
    } else {
        if (knet_close(b->file.remote) != 0)
            err(EX_IOERR, "Error closeing remote BIM file '%s'", b->file_name); 
    }


    free(b->bim_header->md_line_lens);
    free(b->bim_header);
    free(b->gqt_header);

    free(b->file_name);

    b = NULL;
}

struct bim_file_header *read_bim_file_header(struct bim_file *b)
{
    assert(b->type == BIM_REMOTE);

    struct bim_file_header *h = (struct bim_file_header *) 
            malloc(sizeof(struct bim_file_header));

    //if (fseek(b->file, sizeof(struct gqt_file_header), SEEK_SET))
    if (knet_seek(b->file.remote, sizeof(struct gqt_file_header), SEEK_SET))
        err(EX_IOERR,
            "Error seeking to header in BIM file '%s'.",
            b->file_name);

    //size_t fr = fread(&(h->u_size), sizeof(uint64_t), 1, b->file);
    size_t fr = knet_read(b->file.remote, &(h->u_size), 1 * sizeof(uint64_t));
    //check_file_read(b->file_name, b->file, 1, fr);
    check_remote_file_read(b->file_name, 1 * sizeof(uint64_t), fr);

    //fr = fread(&(h->c_size), sizeof(uint64_t), 1, b->file);
    fr = knet_read(b->file.remote, &(h->c_size), 1 * sizeof(uint64_t));
    //check_file_read(b->file_name, b->file, 1, fr);
    check_remote_file_read(b->file_name, 1 * sizeof(uint64_t), fr);

    //fr = fread(&(h->h_size), sizeof(uint64_t), 1, b->file);
    fr = knet_read(b->file.remote, &(h->h_size), 1 * sizeof(uint64_t));
    //check_file_read(b->file_name, b->file, 1, fr);
    check_remote_file_read(b->file_name, 1 * sizeof(uint64_t), fr);

    h->md_line_lens = (uint64_t *) 
            malloc(b->gqt_header->num_variants*sizeof(uint64_t));
    if (!(h->md_line_lens))
        err(EX_OSERR, "malloc error");

    /*
    fr = fread(h->md_line_lens,
               sizeof(uint64_t),
               b->gqt_header->num_variants,
               b->file);
    */
    fr = knet_read(b->file.remote,
                    h->md_line_lens,
                    sizeof(uint64_t) * b->gqt_header->num_variants);

    //check_file_read(b->file_name, b->file, b->gqt_header->num_variants, fr);
    check_remote_file_read(b->file_name,
                           b->gqt_header->num_variants * sizeof(uint64_t),
                           fr);

    return h;
}

struct bim_file *open_bim_file(char *file_name)
{
    struct bim_file *b = (struct bim_file *) malloc(sizeof(struct bim_file));

    b->type = BIM_REMOTE;
    b->file_name = strdup(file_name);
    //b->file = fopen(file_name,"rb+");
    b->file.remote = knet_open(file_name,"rb+");

    if (!(b->file.remote))
        err(EX_NOINPUT, "Cannot open BIM file \"%s\"", file_name);

    b->gqt_header = read_remote_gqt_file_header(b->file_name, b->file.remote);

    if ( !((b->gqt_header->marker[0] == 'G') &&
           (b->gqt_header->marker[1] == 'Q') && 
           (b->gqt_header->marker[2] == 'T')) )
        errx(EX_NOINPUT, "File '%s' is not a GQT file.", file_name);

    if (b->gqt_header->type != 'b')
        errx(EX_NOINPUT, "File '%s' is not a BIM file.", file_name);

    b->bim_header = read_bim_file_header(b);

    //b->data_start = ftell(b->file);
    b->data_start = knet_tell(b->file.remote);

    return b;
}

void update_bim_file_header(uint64_t u_size,
                            uint64_t c_size,
                            uint64_t h_size,
                            struct bim_file *b)
{
    assert(b->type == BIM_LOCAL);

    if (fseek(b->file.local, sizeof(struct gqt_file_header), SEEK_SET))
        err(EX_IOERR,
            "Error seeking to header in BIM file '%s'.",
            b->file_name);

    b->bim_header->u_size = u_size;
    b->bim_header->c_size = c_size;
    b->bim_header->h_size = h_size;

    if (fwrite(&u_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);

    if (fwrite(&c_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);

    if (fwrite(&h_size, sizeof(uint64_t), 1, b->file.local) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);
}

void seek_bim_to_data(struct bim_file *b)
{
    if (b->type == BIM_LOCAL) { 
       if (fseek(b->file.local, b->data_start, SEEK_SET))
            err(EX_IOERR,
                "Error seeking to data in BIM file '%s'.",
                b->file_name);
    } else {
       if (knet_seek(b->file.remote, b->data_start, SEEK_SET) == -1)
            err(EX_IOERR,
                "Error seeking to data in remote BIM file '%s'.",
                b->file_name);
    }
}
