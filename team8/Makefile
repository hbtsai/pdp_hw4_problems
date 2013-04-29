EXEC=a.out
CC=gcc
CXX=g++
DBG_CFLAGS= -D_DEBUG -g 
OBJ=source.o
OMP_CFLAGS=-fopenmp 
CFLAGS= -pthread -fno-builtin -std=c++98  $(DBG_CFLAGS)  $(OMP_CFLAGS)
CXXFLAGS= $(CFLAGS)

all: $(OBJ)
	$(CXX) $(CFLAGS) $^ -o $(EXEC)

%.o: %.c %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

.PHONY: clean

clean:
	rm -f *.o $(EXEC)
