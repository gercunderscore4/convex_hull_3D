CC=g++
CFLAGS=-Wall
LFLAGS=-lm

.PHONY: all clean cleanall

all: ch3d.o main

main: ch3d.o main.cpp
	$(CC) $(CFLAGS) ch3d.o main.cpp -o $@ $(LFLAGS)

ch3d.o: ch3d.h ch3d.cpp
	$(CC) $(CFLAGS) -c ch3d.cpp

clean:
	rm -rf *.o

cleanall:
	rm -rf *.o
	rm -f 'main'
