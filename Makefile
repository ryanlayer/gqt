BIN=bin
OBJ=obj

all: 
	@touch src/parse_q.yy.c
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE)

test: unit_test functional_test

unit_test:
	cd test/unit; $(MAKE); cd ../..

functional_test:
	cd test/func; ./functional_tests.sh

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
