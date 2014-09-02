/*
 * quickFile.h
 *
 *  Created on: Aug 27, 2014
 *      Author: nek3d
 */

#ifndef QUICKFILE_H_
#define QUICKFILE_H_

#include <stddef.h>

/* This will read and store a newline delimited file as an array of char *'s,
 * with a corresponding array of line lengths. */

struct QuickFileInfo {
	char *mainBuf;
	char ** lines;
	size_t *lineLens; /* store length of each line, because strlen is actually quite slow. */
	size_t numLines;
	size_t fileLen;
};

/* Just pass a file name and a pointer to the above struct into this method */
void quickFileInit(char *filename, struct QuickFileInfo *qFile);

/* Don't forget to pass the struct back to this method for memory clean-up */
void quickFileDelete(struct QuickFileInfo *qFile);



#endif /* QUICKFILE_H_ */
