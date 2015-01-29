#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>

int main(int argc, char **argv)
{

    if (argc != 2){
        fprintf(stderr,
                "usage:\t%s <input file>\n",
                argv[0]);
        return 1;
    }

    char *file_name = argv[1];
    FILE *ptr_myfile;
    uint64_t my_record;
    //uint32_t my_record;

    ptr_myfile=fopen(file_name,"rb");
    if (!ptr_myfile) {
        printf("Unable to open file!");
        return 1;
    }

    while ( fread(&my_record,sizeof(uint64_t),1,ptr_myfile) == 1) {
        printf("%" PRIu64 "\n",my_record);
    }

    /*
    while ( fread(&my_record,sizeof(uint32_t),1,ptr_myfile) == 1) {
        printf("%" PRIu32 "\n",my_record);
    }
    */


    fclose(ptr_myfile);
    return 0;
}

