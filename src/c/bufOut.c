/*
 * bufOut.c
 *
 *  Created on: Sep 1, 2014
 *      Author: nek3d
 */

#include "bufOut.h"
#include <string.h>

static int OutputBufferMainBufSize = 1048576; /* 1 megabyte */

/* Allocate and initialize buffer space and members.
 * Pass stdout if not using a file. */
void initOutBuf(struct OutputBuffer *outBuf, FILE *outFile) {
	if ((outBuf->mainBuf = (char *)malloc(OutputBufferMainBufSize)) == NULL) {
		fprintf(stderr, "Error: Failure to allocate output buffer. Exiting...");
		exit(1);
	}
	memset(outBuf->mainBuf, 0, OutputBufferMainBufSize);
	outBuf->outFile = (outFile == NULL ? stdout : outFile);
	outBuf->currPos = 0;
}

/* Send a string and it's length to the output buffer */
void appendOutBuf(struct OutputBuffer *outBuf, char *data, size_t dataLen) {
	if (dataLen + outBuf->currPos >= OutputBufferMainBufSize * 3 / 4) {
		flushOutBuf(outBuf);
	}
	memcpy(outBuf->mainBuf + outBuf->currPos, data, dataLen);
	outBuf->currPos += dataLen;
}


/* deallocate buffer space. Will do one final flush,
 * forced regardless of output size. */
void freeOutBuf(struct OutputBuffer *outBuf) {
	flushOutBuf(outBuf);
	free(outBuf->mainBuf);
}



/* This is just a helper function for the above functions.
 * Flush the contents of the buffer to actual output. */
void flushOutBuf(struct OutputBuffer *outBuf) {
	fprintf(outBuf->outFile, "%s", outBuf->mainBuf);
	memset(outBuf->mainBuf, 0, OutputBufferMainBufSize);
	outBuf->currPos = 0;
}




/* need to calculate the number of digits in the value of R[I] */
void appendIntegerToOutBuf(struct OutputBuffer *outBuf, int data) {


	int copyVal = data;
	int numChars = 0;
	char digits[10];

	if (data < 0) {
		/* Negative number. Include space for minus sign. */
		copyVal = data * -1;
		numChars = 1;
	}
	do {
		numChars++;
		copyVal /= 10;
	}
	while (copyVal > 0);


	if (numChars + outBuf->currPos >= OutputBufferMainBufSize * 3 / 4) {
		flushOutBuf(outBuf);
	}
	sprintf(outBuf->mainBuf + outBuf->currPos, "%u", data);
	outBuf->currPos += numChars;
}
