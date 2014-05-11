BIN=../../bin
OBJ=../../obj


LIB=genotq.o \
    timer.o
OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=and_speed_test \
	get_common_genot_ind \
	get_common_genot_var \
	get_gt_from_ubin \
	int_to_ubin \
	or_gt_fields_from_ubin \
	or_gt_records_from_ubin \
	test_int_genotype_reader \
	or_gt_records_from_uint

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(OBJ)/%.o: %.c %.h
	gcc -c -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	gcc -o $@ $^ 
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
