BIN=../bin
OBJ=../obj
HTS_ROOT=$(HOME)/src/htslib
#HDF5_ROOT=$(HOME)/src/hdf5-1.8.7-mac-intel-x86_64-static
#CFLAGS=-g -D_FILE_OFFSET_BITS=64
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64
#CFLAGS=-Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
#CFLAGS+=-mavx -mavx2
#CFLAGS+=-Dtime_wahbm_pca_by_name
#CFLAGS+=-Dtime_count_range_records_compressed_in_place_wahbm
#CFLAGS+=-Dtime_count_range_records_in_place_wahbm
#CFLAGS+=-Dtime_sum_range_records_in_place_wahbm
#CFLAGS+=-Dtime_convert_file_by_name_bcf_to_wahbm_bim
#CFLAGS+=-Dtime_sort_rotate_gt_md
GRABIX_ROOT=$(HOME)/src/grabix
SQLITE_ROOT=$(HOME)/src/sqlite-amalgamation-3080701
CC=gcc
LEX=flex

LIB= timer.o \
     misc.o \
     query.o \
     parse_q.yy.o \
     parse_q.o \
     ped.o \
     pq.o \
     bcf.o \
     quick_file.o \
     output_buffer.o \
     pthread_pool.o \
     genotq.o \
     convert.o \
     plt.o \
     wah.o \
     wahbm.o \
     wahbm_in_place.o \
     wahbm_compressed_in_place.o \
     ubin.o \
     view.o \
     sort.o \
     count.o \
     sandbox.o

.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=read_binary_uints \
	gqt 


LIST=sqlite parse_q $(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

parse_q:
	$(LEX) -oparse_q.yy.c parse_q.l

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
		-ll -ldl -lz -lm  -pthread
