BIN=../../bin
OBJ=../../obj
CFLAGS=-g


LIB=genotq.o \
    timer.o
OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=and_speed_test \
     	read_binary_uints \
	vcf_to_plt \
	or_records_plt \
	get_common_genot_ind \
	get_common_genot_var \
	get_gt_from_ubin \
	plt_to_ubin \
	ubin_to_wah \
	or_gt_fields_from_ubin \
	or_gt_records_from_ubin \
	test_int_genotype_reader \
	or_gt_records_from_uint \
	or_gt_fields_from_uint \
	sort_plt_by_col

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(OBJ)/%.o: %.c %.h
	gcc -c $(CFLAGS) -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	gcc $(CFLAGS) -o $@ $^ 
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
