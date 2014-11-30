## build for cluster library

COMPILER = gcc 
LIBS =  -lm
CFLAGS = #-std=c89 -pedantic -W -Wall -Wstrict-prototypes -Wunreachable-code  -Wwrite-strings -Wpointer-arith -Wbad-function-cast -Wcast-align -Wcast-qual #-fPIC -g -c -Wall -o
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
	$(COMPILER) -shared -Wl,-soname,libhcluster.so -o libhcluster.so \
	cluster.o   -lc	
