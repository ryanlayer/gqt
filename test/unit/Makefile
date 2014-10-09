#ruby ../auto/generate_test_runner.rb
UNITY_ROOT=./unity
HTS_ROOT=$(HOME)/src/htslib
HDF5_ROOT=$(HOME)/src/hdf5-1.8.7-mac-intel-x86_64-static
CC=gcc
SRC=../../src/c/
TEST_SRC=$(wildcard *.c)
TEST_NAMES=$(TEST_SRC:.c=)
PREFIX=Runner_
RUNNERS=$(addprefix $(PREFIX),$(TESTS))


all: clean $(TEST_NAMES) test

%: %.c
	#ruby $(UNITY_ROOT)/auto/generate_test_runner.rb $@.c
	ruby $(UNITY_ROOT)/generate_test_runner.rb $@.c
	$(CC) -Werror -g -o $@ \
	    -I$(HDF5_ROOT)/include \
	    -L$(HDF5_ROOT)/lib \
	    -lz -lhdf5 -lsz \
	    -I$(UNITY_ROOT) \
	    -I$(SRC) \
	    -I$(HTS_ROOT) \
	    -DTEST $(UNITY_ROOT)/unity.c \
	    -DSAMTOOLS=1 \
	    $(HTS_ROOT)/libhts.a \
	    -mavx -mavx2 -g  -pthread \
	    $(SRC)/genotq.c \
	    $(SRC)/pthread_pool.c \
	    $(SRC)/timer.c \
	    $(SRC)/plt.c \
	    $(SRC)/pq.c \
	    $(SRC)/bcf.c \
	    $(SRC)/h5.c \
	    $(SRC)/ubin.c \
	    $(SRC)/wah.c \
	    $(SRC)/wahbm.c \
	    $(SRC)/wahbm_in_place.c \
	    $(SRC)/wahbm_compressed_in_place.c \
	    $@.c $@_Runner.c -lm 
	rm $@_Runner.c

test:
	$(foreach var,$(TEST_NAMES),./$(var);)

clean:
	rm -f $(TEST_NAMES)
	rm -f *_Runner.c
