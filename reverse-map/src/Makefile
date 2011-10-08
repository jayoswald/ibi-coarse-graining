### BUILD PARAMETERS ##########################################################
CC=g++
EXE=reverse-map
### FIXED PARAMETERES #########################################################
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
CGPOST=../../cg-post/src
###############################################################################

CFLAGS=-c -std=c++0x -I${CGPOST}
LFLAGS=-static

all: CFLAGS+=-O2
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
