#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int int_comp(const void *a,const void *b) 
{
    return ( *(int*)a - *(int*)b );
}

int main(int argc, char **argv)
{

    if (argc < 7) {
        fprintf(stderr,
                "usage:\t%s <indi file> <num ind> <num variants> "
                "<num inds in test> <CSV test inds> <verbose>\n", 
                argv[0]);
        return 1;
    }

    char const* const file_name = argv[1];
    FILE* file = fopen(file_name, "r");
    int num_ind = atoi(argv[2]); 
    int num_var = atoi(argv[3]); 
    int num_test = atoi(argv[4]); 
    char *test_str = argv[5];
    int verbose = atoi(argv[6]);

    int I[num_test];
    char *pch;
    pch = strtok(test_str,",");
    int i;
    for (i = 0; i < num_test; ++i){
        I[i] = atoi(pch);
        pch = strtok(NULL,",");
    }

    qsort(I,num_test,sizeof(int),int_comp) ;

    /* 
    for (i = 0; i < num_test; ++i)
        printf("%d ", I[i]);
    printf("\n");
    */


    char line[2*num_var+1];

    int common_var[num_var];
    for (i = 0; i < num_var; ++i) {
        common_var[i] = 1;
    }

    int var_i,
        line_num = 0,
        curr_I = 0;
    while (fgets(line, sizeof(line), file)) {
        if (line_num == I[curr_I]) {
            var_i = 0;
            pch = strtok(line," ");
            while (pch != NULL) {
                common_var[var_i] = common_var[var_i] & (atoi(pch)>0);
                pch = strtok(NULL," ");
                ++var_i;
            }
            ++curr_I;
        }

        if (curr_I == num_test)
            break;

        ++line_num;

    }
    fclose(file);

    if (verbose > 0) {
        for (i = 0; i < num_var; ++i){
            if (common_var[i] > 0)
                printf("%d\n", i);
        }
    }

    return 0;
}
