#include <stdio.h>
#include <stdlib.h>
#include "genotq.h"

int main(int argc, char **argv)
{
    char *file_name = argv[1];
    char num_gt = atoi(argv[2]);

    init_int_genotype_reader(file_name, num_gt);

    int line_num = 0, gt_num = 0, gt = 0;

    while (get_next_int_genotype(&line_num, &gt_num, &gt) == 0)
        printf("%d\t%d\t%d\n", line_num, gt_num, gt);


    destroy_int_genotype_reader();
}
