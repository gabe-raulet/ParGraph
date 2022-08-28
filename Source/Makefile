CFLAGS=-fopenmp -Wno-unknown-pragmas -g
CC=gcc-11
INCLUDES=-I.
OBJS=spmat.o graph.o minheap.o bitmap.o
PROG=main

D?=0

ifeq ($(D), 1)
CFLAGS+=-O0 -fsanitize=address -fno-omit-frame-pointer -DDEBUG
else
CFLAGS+=-O2
endif

%.o: %.c
	@echo CC -c CFLAGS INCLUDES $<
	@$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ $<

all: $(PROG)

main: main.o $(OBJS)
	@echo CC CFLAGS INCLUDES -o $@ $^
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -rf *.dSYM *.o tests/*.o *.a a.out $(PROG) $(TEST)

spmat.o: spmat.h
graph.o: graph.h spmat.h array.h bitmap.h
minheap.o: minheap.h array.h
bitmap.o: bitmap.h array.h
