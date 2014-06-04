BIN=../../bin
OBJ=../../obj
CFLAGS=-g


LIB=genotq.o \
    timer.o
OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	vcf_to_plt \
	plt_to_ubin \
	ubin_to_wah \
	sort_plt_by_col \
	gt_records_plt \
	gt_records_wah \
	gt_records_ubin

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(OBJ)/%.o: %.c %.h
	gcc -c $(CFLAGS) -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	gcc $(CFLAGS) -o $@ $^ 
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
