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
               "usage:\t%s <vari file> <num ind> <num variants> "
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

    char line[2*num_ind+1];


    int common_var[num_var];
    for (i = 0; i < num_var; ++i) {
        common_var[i] = 1;
    }

    int ind_i,
        var_i = 0,
        curr_I = 0;
    while (fgets(line, sizeof(line), file)) {
        ind_i = 0;
        curr_I = 0;
        pch = strtok(line," ");
        while (pch != NULL) {
            //printf("%d %d %d\n", ind_i, curr_I, I[curr_I]);
            if (ind_i == I[curr_I]) {
                /*
                printf("\tvar_i:%d common_var[var_i]:%d\tpch:%d\tc&p:%d\t", 
                        var_i, 
                        common_var[var_i], 
                        atoi(pch),
                        (common_var[var_i] & atoi(pch))
                        );
                */
                common_var[var_i] = common_var[var_i] & (atoi(pch)>0);
                //printf("%d\n", common_var[var_i]);
                ++curr_I;
            }
            //printf("\n");

            if (curr_I == num_test)
                break;

            pch = strtok(NULL," ");
            ind_i+=1;
        }
        ++var_i;
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
