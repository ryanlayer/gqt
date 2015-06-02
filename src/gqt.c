#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROGRAM_NAME  "gqt"
#define VERSION "0.2.1"


int gqt_help(int argc, char **argv);
int convert(int argc, char **argv);
int view(int argc, char **argv);
int gt(int argc, char **argv);
int sort(int argc, char **argv);
int count(int argc, char **argv);
//int sum(int argc, char **argv);
int query(int argc, char **argv);
int sandbox(int argc, char **argv);
int misc(int argc, char **argv);
int pop(char *op, int argc, char **argv);

int main(int argc, char **argv)
{
    if (argc < 2) return gqt_help(argc, argv);

    char *cmd = argv[1];

    if (strcmp(cmd,"convert") == 0) return convert(argc-2, argv+2);
    else if (strcmp(cmd,"view") == 0) return view(argc-2, argv+2);
    //else if (strcmp(cmd,"gt") == 0) return gt(argc-2, argv+2);
    //else if (strcmp(cmd,"sort") == 0) return sort(argc-2, argv+2);
    //else if (strcmp(cmd,"count") == 0) return count(argc-2, argv+2);
    //else if (strcmp(cmd,"sum") == 0) return sum(argc-2, argv+2);
    else if (strcmp(cmd,"pca-shared") == 0) 
        return pop("pca-shared", argc-1, argv+1);
    else if (strcmp(cmd,"gst") == 0) return pop("gst", argc-1, argv+1);
    else if (strcmp(cmd,"fst") == 0) return pop("fst", argc-1, argv+1);
    else if (strcmp(cmd,"calpha") == 0) return pop("calpha", argc-1, argv+1);
    else if (strcmp(cmd,"query") == 0) return query(argc-1, argv+1);
    else if (strcmp(cmd,"sandbox") == 0) return sandbox(argc-2, argv+2);
    else if (strcmp(cmd,"misc") == 0) return misc(argc-2, argv+2);
    else {
        printf("Unknown command\n");
        return gqt_help(argc, argv);
    }
}

int gqt_help(int argc, char **argv)
{
    printf("%s, v%s\n"
           "usage:   %s <command> [options]\n"
           "         convert    Convert between file types\n"
           "         query      Query the index\n"
           "         pca-shared Compute the similarity matrix for PCA base\n" 
           "                    on the number of shared non-reference loci.\n" 
           "         calpha     Calculate C-alpha paramters (Neal 2011)\n" 
           "         gst        Calculate Gst statistic (Neil 1973)\n"
           "         fst        Calculate Fst statistic " 
                                "(Weir and Cockerham 1984)\n",
            PROGRAM_NAME, VERSION,
            PROGRAM_NAME);
    return 0;
}
