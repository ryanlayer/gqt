#ruby ../auto/generate_test_runner.rb
UNITY_ROOT=/Users/rl6sf/src/Unity-master
CC=gcc
SRC=../c
TEST_SRC=$(wildcard *.c)
TEST_NAMES=$(TEST_SRC:.c=)
PREFIX=Runner_
RUNNERS=$(addprefix $(PREFIX),$(TESTS))


all: clean $(TEST_NAMES) test

%: %.c
	ruby $(UNITY_ROOT)/auto/generate_test_runner.rb $@.c
	$(CC) -Werror -g -o $@ -I$(UNITY_ROOT)/src -I$(SRC) -DTEST $(UNITY_ROOT)/src/unity.c $(SRC)/genotq.c $@.c $@_Runner.c -lm 
	rm $@_Runner.c

test:
	$(foreach var,$(TEST_NAMES),./$(var);)

clean:
	rm -f $(TEST_NAMES)
	rm -f *_Runner.c
