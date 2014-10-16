CC=gcc
CFLAGS=-Wall -Wextra
LFLAGS=-lm

.PHONY: all clean cleanall

all: ch3d.o main

main: ch3d.o main.c
	$(CC) $(CFLAGS) ch3d.o main.c -o $@ $(LFLAGS)

ch3d.o: ch3d.h ch3d.c
	$(CC) $(CFLAGS) -c ch3d.c

clean:
	rm -rf *.o

cleanall:
	rm -rf *.o
	rm -f 'main'
