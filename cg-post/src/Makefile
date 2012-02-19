### BUILD PARAMETERS ##########################################################
CC=g++
EXE=cg-post
### FIXED PARAMETERES #########################################################
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
###############################################################################

CFLAGS=-c -std=c++0x -fopenmp
LFLAGS=-static -fopenmp

all: CFLAGS+=-O3
all: LFLAGS+=-O
all: $(SRC) $(EXE)

debug: CFLAGS+=-g -ggdb
debug: LFLAGS+=-g -ggdb
debug: $(SRC) $(EXE) 

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	@rm -f *.o
