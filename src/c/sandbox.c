/*
 * sandbox.c
 *
 *  Created on: Sep 2, 2014
 *      Author: nek3d
 */

#include <unistd.h>
#include "output_buffer.h"
#include "quick_file.h"
#include "timer.h"


int sandbox(int argc, char **argv)
{

	struct output_buffer outbuf;
	init_out_buf(&outbuf, NULL);
	append_out_buf(&outbuf, "Neil has ", 9);
	append_integer_to_out_buf(&outbuf, 3112);
	append_out_buf(&outbuf, " points.\n", 9);
	append_out_buf(&outbuf, "beer.\n", 6);
	free_out_buf(&outbuf);
}

/*
void spit_back_file(char *inFile);


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
	spit_back_file(in);
	return 0;
}

void spit_back_file(char *inFile) {

	struct quick_file_info qfile;
	struct output_buffer outbuf;
	int i=0;

	start();
	quick_file_init(inFile, &qfile);
	stop();
	fprintf(stderr, "READING AND PARSING THE FILE TOOK %lu MICROSECONDS.\n", report());

	init_out_buf(&outbuf, NULL);


	start();
	for(; i < qfile.num_lines; ++i) {
		append_out_buf(&outbuf, qfile.lines[i], qfile.line_lens[i]);
		append_out_buf(&outbuf, "\n", 1);
	}
	stop();
	fprintf(stderr, "PRINTING THE FILE TOOK %lu MICROSECONDS.\n", report());

	quick_file_delete(&qfile);
	free_out_buf(&outbuf);
}

*/
