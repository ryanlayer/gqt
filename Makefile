BIN=bin
OBJ=obj

all: 
	@touch src/parse_q.yy.c
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE)

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
