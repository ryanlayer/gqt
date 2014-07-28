#ruby ../auto/generate_test_runner.rb
UNITY_ROOT=/Users/rl6sf/src/Unity-master
UNITY_ROOT=./unity
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
	$(CC) -Werror -g -o $@ -I$(UNITY_ROOT) -I$(SRC) -DTEST $(UNITY_ROOT)/unity.c \
	    $(SRC)/genotq.c \
	    $(SRC)/timer.c \
	    $(SRC)/plt.c \
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
