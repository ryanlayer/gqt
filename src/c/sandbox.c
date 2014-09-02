/*
 * sandbox.c
 *
 *  Created on: Sep 2, 2014
 *      Author: nek3d
 */

#include <unistd.h>
#include "bufOut.h"
#include "quickFile.h"

void spitBackFile(char *inFile);

int sandbox(int argc, char **argv)
{
	int c = 0;
	char *in = NULL;
	int minNumArgs = 2;

	if (argc < minNumArgs) {
		fprintf(stderr, "Error: sandbox needs at least %d arguments. Exiting...\n", minNumArgs);
		exit(1);
	}

	in = argv[1];
	spitBackFile(in);
	return 0;
}

void spitBackFile(char *inFile) {

	struct QuickFileInfo qFile;
	struct OutputBuffer outBuf;
	int i=0;

	quickFileInit(inFile, &qFile);

	initOutBuf(&outBuf, NULL);


	for(; i < qFile.numLines; ++i) {
		appendOutBuf(&outBuf, qFile.lines[i], qFile.lineLens[i]);
		appendOutBuf(&outBuf, "\n", 1);
	}
	quickFileDelete(&qFile);
	freeOutBuf(&outBuf);
}
