#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

int sort_help();
int plt_field_freq(char *in, char *out);

struct pair
{
    unsigned int a,b;
};

int compare (const void *a, const void *b)
{
    return ((struct pair *)a)->b - ((struct pair *)b)->b; 
}
int sort(int argc, char **argv)
{
    if (argc < 2) return sort_help();

    int c;
    char *in, *out;
    unsigned int num_fields, num_records;
    int i_is_set = 0, 
        o_is_set = 0; 

    while ((c = getopt (argc, argv, "hi:o:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in = optarg;
                break;
            case 'o':
                o_is_set = 1;
                out = optarg;
                break;
            case 'h':
                sort_help();
                return 1;
            case '?':
                if ( (optopt == 'i') || 
                     (optopt == 'o') )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                sort_help();
                return 1;
        }
    }

    char *type = argv[0];

    if (i_is_set == 0) {
        printf("Input file is not set\n");
        return sort_help();
    } 
    if (o_is_set == 0) {
        printf("Output file is not set\n");
        return sort_help();
    }

    if (strcmp(type, "plt-field-freq") == 0)  return plt_field_freq(in, out);

    return 1;
}

int sort_help()
{
    printf("usage:   gtq sort <type> -i <input file> -o <output file>\n"
           "         plt-field-freq     Sort a plain text file by field "
                                        "frequency\n"
    );

    return 0;
}

int plt_field_freq(char *in, char *out)
{

    struct plt_file pf = init_plt_file(in);

    char *line = NULL;
    size_t len = 0;
    char *pch;
    ssize_t read;
    long line_len = pf.num_fields*2*sizeof(char);
    int i,j;

    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);

    struct pair *S = (struct pair *)
            calloc(pf.num_fields,sizeof(struct pair));

    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        for (j = 0; j < pf.num_fields; ++j) {
            S[j].a = j;
            S[j].b += (unsigned int)line[j*2] - 48;
        }
    }

    qsort(S, pf.num_fields, sizeof(struct pair), compare);


    // jump past the header
    fseek(pf.file, pf.header_offset, SEEK_SET);

    FILE *f = fopen(out, "w");

    fprintf(f, "%d\n%d\n", pf.num_fields, pf.num_records);
    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        for (j = 0; j < pf.num_fields; ++j) {
            if (j != 0)
                fprintf(f, " ");
            fprintf(f, "%c", line[S[j].a*2]);
        }
        fprintf(f, "\n");
    }

    free(line);

    fclose(pf.file);
    fclose(f);

    return 0;
}
