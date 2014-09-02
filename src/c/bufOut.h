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


struct OutputBuffer {
	char *mainBuf;
	size_t currPos;
	FILE *outFile;
};


/* Allocate and initialize buffer space and members.
 * Pass NULL for stdout*/
void initOutBuf(struct OutputBuffer *outBuf, FILE *outFile);

/* Send a string and it's length to the output buffer */
void appendOutBuf(struct OutputBuffer *outBuf, char *data, size_t dataLen);

/* Send an integer to the output buffer */
void appendIntegerToOutBuf(struct OutputBuffer *outBuf, int data);
/* deallocate buffer space. Will do one final flush,
 * forced regardless of output size. */
void freeOutBuf(struct OutputBuffer *outBuf);


/* This is just a helper function for the above functions.
 * Flush the contents of the buffer to actual output. */
void flushOutBuf(struct OutputBuffer *outBuf);

#endif /* _OUPUT_BUFFER_H_ */
