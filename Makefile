BIN=bin
OBJ=obj

all: 
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src/c; make

test:
	cd src/test; make

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
