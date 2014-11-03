BIN=../../bin
OBJ=../../obj
HTS_ROOT=$(HOME)/src/htslib
#HDF5_ROOT=$(HOME)/src/hdf5-1.8.7-mac-intel-x86_64-static
#CFLAGS=-g
CFLAGS=-O3
#CFLAGS=-Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
#CFLAGS+=-mavx -mavx2
#CFLAGS+=-Dtime_count_range_records_compressed_in_place_wahbm
#CFLAGS+=-Dtime_count_range_records_in_place_wahbm
#CFLAGS+=-Dtime_sum_range_records_in_place_wahbm
#CFLAGS+=-Dtime_convert_file_by_name_bcf_to_wahbm_bim
#CFLAGS+=-Dtime_sort_rotate_gt_md
GRABIX_ROOT=$(HOME)/src/grabix
SQLITE_ROOT=$(HOME)/src/sqlite-amalgamation-3080100
CC=gcc

LIB= timer.o \
     ped.o \
     pq.o \
     bcf.o \
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


LIST=sqlite $(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

sqlite: 
	$(CC) -c $(SQLITE_ROOT)/sqlite3.c -o $(SQLITE_ROOT)/sqlite3.o

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< \
		-I$(HTS_ROOT) -I$(SQLITE_ROOT)


$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		-I$(HTS_ROOT) \
		-DSAMTOOLS=1 \
		$(HTS_ROOT)/libhts.a \
		$(SQLITE_ROOT)/sqlite3.o \
		-ldl -lz -lm  -pthread

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
	
