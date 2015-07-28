#ruby ../auto/generate_test_runner.rb
UNITY_ROOT=./unity
HTS_ROOT=$(HOME)/src/htslib
CC=gcc
SRC=../../src/
TEST_SRC=$(wildcard *.c)
TEST_NAMES=$(TEST_SRC:.c=)
PREFIX=Runner_
RUNNERS=$(addprefix $(PREFIX),$(TESTS))
#AVX=-mavx -mavx2

all: clean $(TEST_NAMES) test

%: %.c
	#ruby $(UNITY_ROOT)/auto/generate_test_runner.rb $@.c
	ruby $(UNITY_ROOT)/generate_test_runner.rb $@.c
	lex -o $(SRC)/parse_q.yy.c $(SRC)/parse_q.l
	$(CC) -Werror -g -o $@ \
	    -lz \
	    -I$(UNITY_ROOT) \
	    -I$(SRC) \
	    -I$(HTS_ROOT) \
	    -DTEST $(UNITY_ROOT)/unity.c \
	    -DSAMTOOLS=1 \
	    $(HTS_ROOT)/libhts.a \
	    $(AVX)  -pthread \
	    $(SRC)/genotq.c \
	    $(SRC)/pthread_pool.c \
	    $(SRC)/timer.c \
	    $(SRC)/plt.c \
	    $(SRC)/vid.c \
	    $(SRC)/pq.c \
	    $(SRC)/bcf.c \
	    $(SRC)/ubin.c \
	    $(SRC)/wah.c \
	    $(SRC)/wahbm.c \
	    $(SRC)/wahbm_in_place.c \
	    $(SRC)/wahbm_compressed_in_place.c \
	    $(SRC)/parse_q.c \
	    $(SRC)/parse_q.yy.c \
	    -ll \
	    $@.c $@_Runner.c -lm 
	rm $@_Runner.c

test:
	$(foreach var,$(TEST_NAMES),./$(var);)

clean:
	rm -f $(TEST_NAMES)
	rm -f *_Runner.c
