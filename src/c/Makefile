BIN=../../bin
OBJ=../../obj
CFLAGS=-g 
GRABIX_ROOT=~/src/grabix
CC=gcc

LIB=genotq.o \
    convert.o \
    gt.o \
    view.o \
    sort.o \
    timer.o


OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	gtq \
	gt_records_grabix 

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	# make the "obj/" and "bin/" directories, if they don't exist
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lc
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*

$(BIN)/gt_records_grabix: $(OBJS) gt_records_grabix.cpp $(GRABIX_ROOT)/grabix.cpp
	g++ -c $(GRABIX_ROOT)/grabix.cpp -o $(OBJ)/grabix.o
	gcc -c $(GRABIX_ROOT)/bgzf.c -o $(OBJ)/bgzf.o
	g++ -g gt_records_grabix.cpp \
		genotq.c \
		$(OBJ)/grabix.o \
		$(OBJ)/bgzf.o \
 		-I $(GRABIX_ROOT)  \
		-lstdc++ -lz \
		-o $(BIN)/gt_records_grabix
	
