#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sysexits.h>

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

    b->file = fopen(file_name,"wb");
    if (!b->file)
        err(EX_CANTCREAT, "Cannot create BIM file \"%s\"", file_name);

    if (fwrite(b->gqt_header,
               sizeof(struct gqt_file_header),
               1,
               b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&u_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&c_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(&h_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    if (fwrite(md_line_lens,
               sizeof(uint64_t),
               num_variants,
               b->file) != num_variants)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", file_name);

    b->data_start = ftell(b->file);

    return b;
}

void destroy_bim_file(struct bim_file *b)
{
    if (fclose(b->file))
        err(EX_IOERR, "Error closeing BIM file '%s'", b->file_name); 

    free(b->bim_header->md_line_lens);
    free(b->bim_header);
    free(b->gqt_header);

    free(b->file_name);


    b = NULL;
}

struct bim_file_header *read_bim_file_header(struct bim_file *b)
{
    struct bim_file_header *h = (struct bim_file_header *) 
            malloc(sizeof(struct bim_file_header));

    if (fseek(b->file, sizeof(struct gqt_file_header), SEEK_SET))
        err(EX_IOERR,
            "Error seeking to header in BIM file '%s'.",
            b->file_name);

    size_t fr = fread(&(h->u_size), sizeof(uint64_t), 1, b->file);
    check_file_read(b->file_name, b->file, 1, fr);

    fr = fread(&(h->c_size), sizeof(uint64_t), 1, b->file);
    check_file_read(b->file_name, b->file, 1, fr);

    fr = fread(&(h->h_size), sizeof(uint64_t), 1, b->file);
    check_file_read(b->file_name, b->file, 1, fr);

    h->md_line_lens = (uint64_t *) 
            malloc(b->gqt_header->num_variants*sizeof(uint64_t));
    if (!(h->md_line_lens))
        err(EX_OSERR, "malloc error");

    fr = fread(h->md_line_lens,
               sizeof(uint64_t),
               b->gqt_header->num_variants,
               b->file);
    check_file_read(b->file_name, b->file, b->gqt_header->num_variants, fr);

    return h;
}

struct bim_file *open_bim_file(char *file_name)
{
    struct bim_file *b = (struct bim_file *) malloc(sizeof(struct bim_file));

    b->file_name = strdup(file_name);
    b->file = fopen(file_name,"rb+");

    if (!(b->file))
        err(EX_NOINPUT, "Cannot open BIM file \"%s\"", file_name);

    b->gqt_header = read_gqt_file_header(b->file_name, b->file);

    if ( !((b->gqt_header->marker[0] == 'G') &&
           (b->gqt_header->marker[1] == 'Q') && 
           (b->gqt_header->marker[2] == 'T')) )
        errx(EX_NOINPUT, "File '%s' is not a GQT file.", file_name);

    if (b->gqt_header->type != 'b')
        errx(EX_NOINPUT, "File '%s' is not a BIM file.", file_name);

    b->bim_header = read_bim_file_header(b);

    b->data_start = ftell(b->file);

    return b;
}

void update_bim_file_header(uint64_t u_size,
                            uint64_t c_size,
                            uint64_t h_size,
                            struct bim_file *b)
{
    if (fseek(b->file, sizeof(struct gqt_file_header), SEEK_SET))
        err(EX_IOERR,
            "Error seeking to header in BIM file '%s'.",
            b->file_name);

    b->bim_header->u_size = u_size;
    b->bim_header->c_size = c_size;
    b->bim_header->h_size = h_size;

    if (fwrite(&u_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);

    if (fwrite(&c_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);

    if (fwrite(&h_size, sizeof(uint64_t), 1, b->file) != 1)
        err(EX_IOERR, "Error writing header to BIM file \"%s\"", b->file_name);
}

void seek_bim_to_data(struct bim_file *b)
{
   if (fseek(b->file, b->data_start, SEEK_SET))
        err(EX_IOERR,
            "Error seeking to data in BIM file '%s'.",
            b->file_name);
}
