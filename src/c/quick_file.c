/*
 * quickFile.c
 *
 *  Created on: Aug 27, 2014
 *      Author: nek3d
 */

#include "quick_file.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "genotq.h"


void quick_file_init(char *filename, struct quick_file_info *qfile) {
    FILE *fp = NULL;
    uint64_t i = 0;
    uint64_t pos = 0;

    /* Open the file */
    if ((fp = fopen(filename, "rb")) == NULL) {
        fprintf(stderr,
                "Error: Unable to open file %s. Exiting....\n",
                filename);
        exit(1);
    }

    /* 
     * The file is :
     * uncompressed size  ( sizeof(size_t))
     * compressed size    ( sizeof(size_t))
     * header size        ( sizeof(size_t))
     * compressed data 
     */
    uint64_t u_size, c_size, h_size;
    size_t s = fread(&u_size, sizeof(uint64_t), 1, fp);
    s = fread(&c_size, sizeof(uint64_t), 1, fp);
    s = fread(&h_size, sizeof(uint64_t), 1, fp);

    //fprintf(stderr, "u_size:%u\tc_size:%u\th_size:%u\n",
            //u_size, c_size, h_size);


    /* allocate inflate state */
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    int ret = inflateInit(&strm);
    if (ret != Z_OK) {
        fprintf(stderr, "error: Cannot init stream\n");
        exit(1);
    }

    // in_buf will hold compressed data and qfile->main_buf uncompressed
    //unsigned char *in_buf = (unsigned char *)
            //malloc(sizeof(unsigned char) * CHUNK);
    unsigned char in_buf[CHUNK];
            

    char *final_out_buf = (char *) malloc(u_size);
    qfile->main_buf = final_out_buf;

    strm.avail_out = u_size;
    strm.next_out = (Bytef *)final_out_buf;

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = fread(in_buf, 1, CHUNK, fp);
        if (ferror(fp)) {
            fprintf(stderr, "error: Cannot read compressed file.\n");
            (void)inflateEnd(&strm);
            exit(1);
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in_buf;

        ret = inflate(&strm, Z_NO_FLUSH);
        assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                exit(1);
        }

    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);

#if 0
    Bytef *c_buf = (Bytef *)malloc(c_len);
    s = fread(c_buf, c_len, 1, fp);
    fclose(fp);

    /* allocate the main buffer to hold the entire file in one
     * contiguous block. */

    int r = uncompress((Bytef *)(qfile->main_buf), &u_len, c_buf, c_len);
    if (r == Z_BUF_ERROR) {
        fprintf(stderr,
                "The buffer dest was not large enough to hold the "
                "uncompressed data.");
        exit(1);
    } else if (r == Z_MEM_ERROR) {
        fprintf(stderr, "Insufficient memory");
        exit(1);
    } else if (r == Z_DATA_ERROR) {
        fprintf(stderr,
                "The compressed data (referenced by source) was corrupted.");
        exit(1);
    }


    free(c_buf);
#endif

    qfile->file_len = u_size;
    qfile->header_len = h_size;

    /* 
     * Count how many lines you have.  Starting after the header
     * 
     */
    qfile->num_lines = 0;
    for(i = qfile->header_len; i < qfile->file_len; i++) {
        if (qfile->main_buf[i] == '\n') qfile->num_lines++;
    }

    /* Special: if the last char in the file is not a newLine,
     * then there was one more line we missed. */
    if (qfile->main_buf[i-1] != '\n') {
        qfile->num_lines++;
    }

    /* Side note: the reason we made a pass just to count newlines,
     * and then are passing through again below, is so that we know
     * how many elements to allocate for an array of newlines,
     * rather than try to write a dynamic array, which has to do
     * reallocations and copies, because that's way slower than just
     * doing a 2nd pass through the data.
     */


    /* allocate an array of char **'s, one for each line. */
    qfile->lines = (char **)malloc(qfile->num_lines *sizeof(char *));
    memset(qfile->lines, 0, qfile->num_lines * sizeof(char *));

    /* also allocate the array of line lengths. */
    qfile->line_lens = (uint64_t *)malloc(qfile->num_lines * sizeof(uint64_t));

    /* Loop through the mainbuf again, turning all newlines into
     * null chars. Put a pointer to the next char on the lines array.
     * Store the calculated length of the line as well. */
    i = 0;
    qfile->lines[0] = qfile->main_buf + qfile->header_len;

    uint64_t prevPos = qfile->header_len;

    for (pos = qfile->header_len; pos < qfile->file_len; pos++) {
        if (qfile->main_buf[pos] == '\n') {
            i++;
            qfile->line_lens[i-1] = pos - prevPos;
            prevPos = pos +1;
            qfile->main_buf[pos] = '\0';

            /* If the newline wasn't the last char in the file,
             * make a new line. */
            if (pos < qfile->file_len-1) {
                    qfile->lines[i] = qfile->main_buf + pos +1;
            }
        }
    }

    /* Special: If the file didn't end with a newline, there
     * is a partial line we missed. */
    if (prevPos < qfile->file_len) {
        qfile->line_lens[i] = qfile->file_len - prevPos;
    }
    /* And we're done! */
}

void quick_file_delete(struct quick_file_info *qfile) {
    /* This function is just for memory cleanup. */
    int pos = 0;

    /* free the array of lines. */
    free(qfile->lines);

    /* Turn all the null chars in the main_buf, except the last null char,
     * back into newlines, so that the main_buf can be freed all at once. */
    for (pos = 0; pos < qfile->file_len -1; pos++) {
        if (qfile->main_buf[pos] == '\0') {
                qfile->main_buf[pos] = '\n';
        }
    }


    /* now we can free the main_buf. */
    free(qfile->main_buf);
}
