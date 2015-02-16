## build for cluster library

COMPILER = gcc 
LIBS =  -lm
CFLAGS = #-std=c89 -pedantic -W -Wall -Wstrict-prototypes -Wunreachable-code  -Wwrite-strings -Wpointer-arith -Wbad-function-cast -Wcast-align -Wcast-qual #
COPTS =  -Wextra -ansi -pedantic -Wmissing-prototypes -Wstrict-prototypes \
    -Wold-style-definition -fPIC -g -c -Wall -fwrapv -O3 -fno-strict-aliasing -march=nocona -O3 -mmmx -msse -pthread -o 

all: cluster

clean:
	rm *.o *.so

dist:
	tar -cf ../cluster.tar *.c *.h Makefile README.md LICENSE 

hclust.o: hclust.c
	$(COMPILER) $(COPTS) hclust.o hclust.c  $(LIBS) $(CFLAGS)

cluster: hclust.o
	$(COMPILER) -shared -Wstrict-prototypes -fno-strict-aliasing -Wl,-soname,libhcluster.so -o libhcluster.so \
	 hclust.o -lc -lm

