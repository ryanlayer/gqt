#include <stdio.h>
#include <stdlib.h>
#include "timer.h"


int main(int argc, char **argv)
{

    int I = atoi(argv[1]);

    int i;
    int *v = (int *) malloc (I*sizeof(int));

    for (i = 1; i < I-1; ++i)
        v[i] = v[i-1] + v[i]; 

    for (i = 1; i < I-1; ++i)
        v[i] = v[i-1] & v[i]; 

    start();
    for (i = 1; i < I-1; ++i)
        v[i] = v[i-1] + v[i]; 
    stop();

    printf("%lu\n", report());

    start();
    for (i = 1; i < I-1; ++i)
        v[i] = v[i-1] & v[i]; 
    stop();

    printf("%lu\n", report());

}
