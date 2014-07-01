#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROGRAM_NAME  "gtq"
#define VERSION "0.0.1"


int gtq_help(int argc, char **argv);
int convert(int argc, char **argv);
int view(int argc, char **argv);

int main(int argc, char **argv)
{
    if (argc < 2) return gtq_help(argc, argv);

    char *cmd = argv[1];

    if (strcmp(cmd,"convert") == 0) return convert(argc-2, argv+2);
    else if (strcmp(cmd,"view") == 0) return view(argc-2, argv+2);
    else {
        printf("Unknown command\n");
        return gtq_help(argc, argv);
    }
}

int gtq_help(int argc, char **argv)
{
    printf("%s, v%s\n"
           "usage:   %s <command> [options]\n"
           "         convert   Convert between file types\n"
           "         view      Display files contents\n",
            PROGRAM_NAME, VERSION,
            PROGRAM_NAME);
    return 0;
}
