## build for cluster library

COMPILER = gcc 
LIBS =  -lm
CFLAGS = #-fPIC -g -c -Wall -o
COPTS =  -fPIC -g -c -Wall -lm -march=nocona -O3 -mmmx -msse -pthread -o 
LIBPATH = 



all: cluster

clean:
	rm *.o *.so

dist:
	tar -cf ../cluster.tar *.c *.h Makefile README.md LICENSE 

cluster.o: cluster.c
	$(COMPILER) $(COPTS) cluster.o cluster.c  $(LIBS) $(CFLAGS)

cluster: cluster.o 
	$(COMPILER) -shared -Wl,-soname,cluster.so -o cluster.so \
	cluster.o   -lc	
