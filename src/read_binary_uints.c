#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <err.h>
#include <sysexits.h>

int main(int argc, char **argv)
{

    if (argc != 3){
        fprintf(stderr,
                "usage:\t%s <size> <input file>\n",
                argv[0]);
        return 1;
    }

    int size = atoi(argv[1]);
    char *file_name = argv[2];
    FILE *ptr_myfile;
    ptr_myfile=fopen(file_name,"rb");
    if (!ptr_myfile)
        err(EX_NOINPUT, "Cannot open \"%s\"", file_name);

    if (!ptr_myfile) {
        printf("Unable to open file!");
        return 1;
    }

    if (size == 32) {
        uint32_t my_record;
        while ( fread(&my_record,sizeof(uint32_t),1,ptr_myfile) == 1) {
            printf("%" PRIu32 "\n",my_record);
        }
    } else if (size == 64) {
        uint64_t my_record;
        while ( fread(&my_record,sizeof(uint64_t),1,ptr_myfile) == 1) {
            printf("%" PRIu64 "\n",my_record);
        }
    }
    fclose(ptr_myfile);
    return 0;
}

