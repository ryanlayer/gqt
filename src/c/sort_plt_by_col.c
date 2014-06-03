#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "genotq.h"

struct pair
{
    unsigned int a,b;
};

int compare (const void *a, const void *b)
{
    return ((struct pair *)a)->b - ((struct pair *)b)->b; 
}

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\%s <options>\n"
        "\t\t-i\tInput variant or individual plain text file name\n"
        , prog
    );
}

int main(int argc, char **argv)
{
    char *prog = argv[0];
    int c;
    char *in_file_name;
    int i_is_set = 0;

    while ((c = getopt (argc, argv, "hi:")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                in_file_name = optarg;
                break;
            case 'h':
                usage(prog);
                return 1;
            case '?':
                if ( optopt == 'i' )
                    fprintf (stderr, "Option -%c requires an argument.\n",
                            optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                usage(prog);
                return 1;
        }
    }

    if (  i_is_set == 0 ) {
        usage(prog);
        return 1;
    }

    struct plt_file pf = init_plt_file(in_file_name);

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

    fseek(pf.file, pf.header_offset, SEEK_SET);

    printf("%d\n%d\n", pf.num_fields, pf.num_records);
    for (i = 0; i < pf.num_records; ++i) {
        read = getline(&line, &len, pf.file);
        for (j = 0; j < pf.num_fields; ++j) {
            if (j != 0)
                printf(" ");
            printf("%c", line[S[j].a*2]);
        }
        printf("\n");
    }

    free(line);

    fclose(pf.file);

    return 0;
}
