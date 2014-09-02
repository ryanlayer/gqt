/*
 * quickFile.c
 *
 *  Created on: Aug 27, 2014
 *      Author: nek3d
 */

#include "quickFile.h"
#include <stdlib.h>
#include <stdio.h>
#include "sys/stat.h"


void quickFileInit(char *filename, struct QuickFileInfo *qFile) {
	struct stat infileStat;
	FILE *fp = NULL;
	size_t i = 0;
	size_t pos = 0;
	size_t prevPos = 0;

	/* Open the file */
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Error: Unable to open file %s. Exiting....\n", filename);
		exit(1);
	}

	/*figure out what size the file is */
	stat(filename, &infileStat);
	qFile->fileLen = infileStat.st_size;


	/* allocate the main buffer to hold the entire file in one
	 * contiguous block. */
	qFile->mainBuf = (char *)malloc(qFile->fileLen +1);
	memset(qFile->mainBuf, 0, qFile->fileLen +1);

	/* read the file into that buffer in one gulp, then close the file */
	if (fread(qFile->mainBuf, 1, qFile->fileLen, fp) < qFile->fileLen) {
		fprintf(stderr, "Error: Unable to read in all of file %s. Exiting...\n ", filename);
		exit(1);
	}
	fclose(fp);



	/* Count how many lines you have. There may be a better way
	 * to do this, but I don't know what it is.  */
	qFile->numLines = 0;
	for(i = 0; i < qFile->fileLen; i++) {
		if (qFile->mainBuf[i] == '\n') qFile->numLines++;
	}
	/* Special: if the last char in the file is not a newLine,
	 * then there was one more line we missed. */
	if (qFile->mainBuf[i-1] != '\n') {
		qFile->numLines++;
	}


	/* Side note: the reason we made a pass just to count newlines,
	 * and then are passing through again below, is so that we know
	 * how many elements to allocate for an array of newlines,
	 * rather than try to write a dynamic array, which has to do
	 * reallocations and copies, because that's way slower than just
	 * doing a 2nd pass through the data.
	 */


	/* allocate an array of char **'s, one for each line. */
	qFile->lines = (char **)malloc(qFile->numLines *sizeof(char *));
	memset(qFile->lines, 0, qFile->numLines * sizeof(char *));

	/* also allocate the array of line lengths. */
	qFile->lineLens = (size_t *)malloc(qFile->numLines * sizeof(size_t));

	/* Loop through the mainbuf again, turning all newlines into
	 * null chars. Put a pointer to the next char on the lines array.
	 * Store the calculated length of the line as well. */
	i = 0;
	qFile->lines[0] = qFile->mainBuf;
	for (pos = 0; pos < qFile->fileLen; pos++) {
		if (qFile->mainBuf[pos] == '\n') {
			i++;
			qFile->lineLens[i-1] = pos - prevPos;
			prevPos = pos +1;
			qFile->mainBuf[pos] = '\0';

			/* If the newline wasn't the last char in the file,
			 * make a new line. */
			if (pos < qFile->fileLen-1) {
				qFile->lines[i] = qFile->mainBuf + pos +1;
			}
		}
	}

	/* Special: If the file didn't end with a newline, there
	 * is a partial line we missed. */
	if (prevPos < qFile->fileLen) {
		qFile->lineLens[i] = qFile->fileLen - prevPos;
	}
	/* And we're done! */
}

void quickFileDelete(struct QuickFileInfo *qFile) {
	/* This function is just for memory cleanup. */
	int pos = 0;

	/* free the array of lines. */
	free(qFile->lines);

	/* Turn all the null chars in the mainBuf, except the last null char,
	 * back into newlines, so that the mainBuf can be freed all at once. */
	for (pos = 0; pos < qFile->fileLen -1; pos++) {
		if (qFile->mainBuf[pos] == '\0') {
			qFile->mainBuf[pos] = '\n';
		}
	}

	/* now we can free the mainBuf. */
	free(qFile->mainBuf);
}





