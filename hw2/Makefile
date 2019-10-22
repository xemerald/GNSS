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

all: hw2_1 hw2_3

hw2_1: $(LIBS) $(LLIBS)
	$(CC) $(CFLAGS) -o ./hw2_1 hw2_1.c $(LIBS) $(LLIBS)

hw2_3: hw2_3.o $(LIBS) $(LLIBS)
	$(FCOMPL) -o hw2_3 hw2_3.o $(LLIBS) $(CPGPLOT_LIB) $(LIBS)
	export PGPLOT_DIR="../data/"
	
hw2_3_matrix: hw2_3_matrix.o $(LIBS) $(LLIBS)
	$(FCOMPL) -o hw2_3_matrix hw2_3_matrix.o $(LLIBS) $(CPGPLOT_LIB) $(LIBS)
	export PGPLOT_DIR="../data/"
	
# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f ./hw2_1 ./hw2_3
