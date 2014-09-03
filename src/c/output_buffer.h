/*
 * bufOut.h
 *
 *  Created on: Sep 1, 2014
 *      Author: nek3d
 */

#ifndef _OUPUT_BUFFER_H_
#define _OUPUT_BUFFER_H_

#include <stdio.h>
#include <stdlib.h>


struct output_buffer {
	char *main_buf;
	size_t curr_pos;
	FILE *out_file;
};


/* Allocate and initialize buffer space and members.
 * Pass NULL for stdout*/
void init_out_buf(struct output_buffer *out_buf, FILE *out_file);

/* Send a string and it's length to the output buffer */
void append_out_buf(struct output_buffer *out_buf, char *data, size_t data_len);

/* Send an integer to the output buffer */
void append_integer_to_out_buf(struct output_buffer *out_buf, int data);
/* deallocate buffer space. Will do one final flush,
 * forced regardless of output size. */
void free_out_buf(struct output_buffer *out_buf);


/* This is just a helper function for the above functions.
 * Flush the contents of the buffer to actual output. */
void flush_out_buf(struct output_buffer *out_buf);

#endif /* _OUPUT_BUFFER_H_ */
