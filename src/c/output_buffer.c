/*
 * bufOut.c
 *
 *  Created on: Sep 1, 2014
 *      Author: nek3d
 */

#include "output_buffer.h"
#include <string.h>

static int output_buffer_main_buf_size = 1048576; /* 1 megabyte */

/* Allocate and initialize buffer space and members.
 * Pass stdout if not using a file. */
void init_out_buf(struct output_buffer *out_buf, FILE *out_file) {
	if ((out_buf->main_buf = (char *)malloc(output_buffer_main_buf_size)) == NULL) {
		fprintf(stderr, "Error: Failure to allocate output buffer. Exiting...");
		exit(1);
	}
	memset(out_buf->main_buf, 0, output_buffer_main_buf_size);
	out_buf->out_file = (out_file == NULL ? stdout : out_file);
	out_buf->curr_pos = 0;
}

/* Send a string and it's length to the output buffer */
void append_out_buf(struct output_buffer *out_buf, char *data, size_t data_len) {
	if (data_len + out_buf->curr_pos >= output_buffer_main_buf_size) {
		flush_out_buf(out_buf);
	}
	memcpy(out_buf->main_buf + out_buf->curr_pos, data, data_len);
	out_buf->curr_pos += data_len;
}


/* deallocate buffer space. Will do one final flush,
 * forced regardless of output size. */
void free_out_buf(struct output_buffer *out_buf) {
	flush_out_buf(out_buf);
	free(out_buf->main_buf);
}



/* This is just a helper function for the above functions.
 * Flush the contents of the buffer to actual output. */
void flush_out_buf(struct output_buffer *out_buf) {
	fprintf(out_buf->out_file, "%s", out_buf->main_buf);
	memset(out_buf->main_buf, 0, output_buffer_main_buf_size);
	out_buf->curr_pos = 0;
}




/* need to calculate the number of digits in the value of R[I] */
void append_integer_to_out_buf(struct output_buffer *out_buf, int data) {


	int copy_val = data;
	int num_chars = 0;
	char digits[10];

	if (data < 0) {
		/* Negative number. Include space for minus sign. */
		copy_val = data * -1;
		num_chars = 1;
	}
	do {
		num_chars++;
		copy_val /= 10;
	}
	while (copy_val > 0);


	if (num_chars + out_buf->curr_pos >= output_buffer_main_buf_size * 3 / 4) {
		flush_out_buf(out_buf);
	}
	sprintf(out_buf->main_buf + out_buf->curr_pos, "%u", data);
	out_buf->curr_pos += num_chars;
}
