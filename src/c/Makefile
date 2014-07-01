BIN=../../bin
OBJ=../../obj
CFLAGS=-g
GRABIX_ROOT=~/src/grabix


LIB=genotq.o \
    convert.o \
    view.o \
    timer.o


OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	vcf_to_plt \
	plt_to_ubin \
	ubin_to_wah \
	sort_plt_by_col \
	gt_records_plt \
	gt_records_wah \
	gt_records_ubin \
	compress_ubin_with_wah \
	gtq \
	gt_records_grabix 

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(OBJ)/%.o: %.c 
	gcc -c $(CFLAGS) -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	gcc $(CFLAGS) -o $@ $^ -lc
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*

$(BIN)/gt_records_grabix: gt_records_grabix.cpp $(GRABIX_ROOT)/grabix.cpp
	g++ -g gt_records_grabix.cpp \
	    	genotq.c \
		$(GRABIX_ROOT)/grabix.cpp \
		$(GRABIX_ROOT)/bgzf.c \
 		-I $(GRABIX_ROOT)  \
		-lstdc++ -lz \
		-o $(BIN)/gt_records_grabix
	
