#
#   THIS FILE IS UNDER RCS - DO NOT MODIFY UNLESS YOU HAVE
#   CHECKED IT OUT USING THE COMMAND CHECKOUT.
#
#

CFLAGS = $(GLOBALFLAGS) -O3 -g -I../include -I/home/benyang/ew/pgplot -flto
LIBS = -lm -lpng -lz -lX11
FCOMPL = gfortran

L = ../lib

LLIBS = $L/xyz2geoc.o $L/matrix.o

CPGPLOT_LIB = /home/benyang/ew/pgplot/libcpgplot.a /home/benyang/ew/pgplot/libpgplot.a

all: hw3_1

hw3_1: hw3_1.o $(LIBS) $(LLIBS)
	$(FCOMPL) -o hw3_1 hw3_1.o $(LLIBS) $(CPGPLOT_LIB) $(LIBS)
	export PGPLOT_DIR="../data/"

hw3_2: hw3_2.o $(LIBS) $(LLIBS)
	$(CC) $(CFLAGS) -o hw3_2 hw3_2.o $(LLIBS) $(LIBS)

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f ./hw2_1 ./hw2_3
