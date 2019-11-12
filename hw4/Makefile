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

all: hw4_2

hw4_2: hw4_2.o $(LIBS) $(LLIBS)
	$(CC) $(CFLAGS) -o hw4_2 hw4_2.o $(LLIBS) $(LIBS)

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f ./hw4_2
