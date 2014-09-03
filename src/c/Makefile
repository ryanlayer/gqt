BIN=../../bin
OBJ=../../obj
#CFLAGS=-g
CFLAGS=-O3
CFLAGS+=-mavx -mavx2
#CFLAGS+=-Dtime_count_range_records_compressed_in_place_wahbm
#CFLAGS+=-Dtime_count_range_records_in_place_wahbm
#CFLAGS+=-Dtime_sum_range_records_in_place_wahbm
#CFLAGS+=-Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
GRABIX_ROOT=~/src/grabix
CC=gcc

LIB= timer.o \
	quickFile.o \
    pthread_pool.o \
    genotq.o \
    convert.o \
    gt.o \
    plt.o \
    wah.o \
    wahbm.o \
    wahbm_in_place.o \
    wahbm_compressed_in_place.o \
    ubin.o \
    view.o \
    sort.o \
    count.o \
    sum.o 


OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	gtq 


LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	# make the "obj/" and "bin/" directories, if they don't exist
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< 

$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lc -lm  -pthread
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*

$(BIN)/gt_records_grabix: $(OBJS) gt_records_grabix.cpp $(GRABIX_ROOT)/grabix.cpp
	g++ -c $(GRABIX_ROOT)/grabix.cpp -o $(OBJ)/grabix.o
	gcc -c $(GRABIX_ROOT)/bgzf.c -o $(OBJ)/bgzf.o
	g++ -g gt_records_grabix.cpp \
		quickFile.c \
		genotq.c \
		wah.c \
		wahbm.c \
		wahbm_in_place.c \
		wahbm_compressed_in_place.c \
		$(OBJ)/grabix.o \
		$(OBJ)/bgzf.o \
 		-I $(GRABIX_ROOT)  \
		-lstdc++ -lz \
		-o $(BIN)/gt_records_grabix
	
