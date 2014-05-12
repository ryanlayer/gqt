#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

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
    unsigned int my_record;

    ptr_myfile=fopen(file_name,"rb");
    if (!ptr_myfile) {
        printf("Unable to open file!");
        return 1;
    }

    while ( fread(&my_record,sizeof(unsigned int),1,ptr_myfile) == 1) {
        printf("%u\n",my_record);
    }

    fclose(ptr_myfile);
    return 0;
}

