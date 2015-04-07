#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

void usage(char *prog)
{
    fprintf(stderr,
        "usage:\t%s<options>\n"
        "\t\t-u\tInput variant or individual uncompressed binary file name\n"
        "\t\t-n\tNumber of fields to OR\n"
        "\t\t-f\tCSV field IDs\n", prog
    );
}


int main(int argc, char **argv)
{

    //{{{ option parsing
    char *prog = argv[0];
    int num_fields,
        num_records,
        r_is_set = 0,
        f_is_set = 0; 
    int c;

    while ((c = getopt (argc, argv, "f:r:")) != -1) {
        switch (c) {
            case 'f':
                f_is_set = 1;
                num_fields = atoi(optarg);
                break;
            case 'r':
                r_is_set = 1;
                num_records = atoi(optarg);
                break;
            case 'h':
                usage(prog);
                return 1;
            case '?':
                if ((optopt == 'f') || (optopt == 'r'))
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

    if (  (f_is_set == 0) || 
          (r_is_set == 0) ) {
        usage(prog);
        exit(EXIT_SUCCESS);
    }
    //}}}

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *pch;
    int gt;

    printf("%d\n%d\n", num_fields, num_records);

    while ((read = getline(&line, &len, stdin)) != -1) {
        if (line[0] != '#') {

            // Skip the first 9 fields
            
            int i;
            pch = strtok(line, "\t");
            for  (i =0; i < 8; ++i) {
                pch = strtok(NULL, "\t");
            }
            
            // get the first genotype
            pch = strtok(NULL, "\t");
            i = 0;
            while (pch != NULL) {
                if ((pch[0] == '0') && (pch[2] == '0'))
                    gt = 0;
                else if ((pch[0] == '1') && (pch[2] == '0'))
                    gt = 1;
                else if ((pch[0] == '0') && (pch[2] == '1'))
                    gt = 1;
                else if ((pch[0] == '1') && (pch[2] == '1'))
                    gt = 2;
                else
                    gt = 3;

                if (i != 0)
                    printf(" ");

                printf("%d", gt);
                pch = strtok(NULL, "\t");
                i+=1;
            }
            printf("\n");
        }
    }

    free(line);
    exit(EXIT_SUCCESS);

}
