BIN=../../bin
OBJ=../../obj
HTS_ROOT=$(HOME)/src/htslib
HDF5_ROOT=$(HOME)/src/hdf5-1.8.7-mac-intel-x86_64-static
#CFLAGS=-g
CFLAGS=-O3
#CFLAGS=-Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
#CFLAGS+=-mavx -mavx2
#CFLAGS+=-Dtime_count_range_records_compressed_in_place_wahbm
#CFLAGS+=-Dtime_count_range_records_in_place_wahbm
#CFLAGS+=-Dtime_sum_range_records_in_place_wahbm
GRABIX_ROOT=~/src/grabix
CC=gcc

LIB= timer.o \
     pq.o \
     bcf.o \
     h5.o \
     quick_file.o \
     output_buffer.o \
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
     sum.o \
     sandbox.o


OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	gqt 


LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< \
		-I$(HDF5_ROOT)/include \
		-I$(HTS_ROOT) 


$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		-I$(HDF5_ROOT)/include \
		-I$(HTS_ROOT) \
		-L$(HDF5_ROOT)/lib \
		-DSAMTOOLS=1 \
		$(HTS_ROOT)/libhts.a \
		-lz -lm  -pthread -lhdf5 -lsz 

$(BIN)/gt_records_grabix: $(OBJS) gt_records_grabix.cpp $(GRABIX_ROOT)/grabix.cpp
	g++ -c $(GRABIX_ROOT)/grabix.cpp -o $(OBJ)/grabix.o
	gcc -c $(GRABIX_ROOT)/bgzf.c -o $(OBJ)/bgzf.o
	g++ -g gt_records_grabix.cpp \
		quick_file.c \
		buf_out.c \
		genotq.c \
		wah.c \
		wahbm.c \
		wahbm_in_place.c \
		wahbm_compressed_in_place.c \
		sandbox.c \
		$(OBJ)/grabix.o \
		$(OBJ)/bgzf.o \
 		-I $(GRABIX_ROOT)  \
		-lstdc++ -lz \
		-o $(BIN)/gt_records_grabix
	
