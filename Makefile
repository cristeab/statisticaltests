# Needs GNU Scientific Library (gsl) installed
# Optimizations are made for Athlon64 processor, change -march switch
#for other processor type or delete this switch

CC=g++
CFLAGS=-Wall -O3 -pipe -march=athlon64 -lgsl
DFLAGS=-ggdb -lgsl

all : nist

debug : nist_d

nist : nist.cpp Statistical_tests.cpp
	$(CC) $(CFLAGS) nist.cpp -o nist

nist_d : nist.cpp Statistical_tests.cpp
	$(CC) $(DFLAGS) nist.cpp -o nist_d

help :
	doxygen Doxyfile

clean :
	rm -f nist nist_d
	rm -f *.orig