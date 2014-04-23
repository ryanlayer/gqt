
all: get_common_genot_ind \
	get_common_genot_var \
	test_int_genotype_reader \
	int_to_ubin \
	get_gt_from_ubin \
	or_gt_records_from_ubin \
	or_gt_fields_from_ubin


get_common_genot_ind: get_common_genot_ind.c
	gcc -o get_common_genot_ind get_common_genot_ind.c
	
get_common_genot_var: get_common_genot_var.c
	gcc -o get_common_genot_var get_common_genot_var.c

test_int_genotype_reader: test_int_genotype_reader.c genotq.c
	gcc -g -o test_int_genotype_reader test_int_genotype_reader.c genotq.c

int_to_ubin: int_to_ubin.c genotq.c
	gcc -g -o int_to_ubin int_to_ubin.c genotq.c

get_gt_from_ubin: get_gt_from_ubin.c genotq.c
	gcc -g -o get_gt_from_ubin get_gt_from_ubin.c genotq.c

or_gt_records_from_ubin: or_gt_records_from_ubin.c genotq.c
	gcc -g -o or_gt_records_from_ubin or_gt_records_from_ubin.c genotq.c

or_gt_fields_from_ubin: or_gt_fields_from_ubin.c genotq.c
	gcc -g -o or_gt_fields_from_ubin or_gt_fields_from_ubin.c genotq.c

