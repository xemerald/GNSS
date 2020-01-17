#
#   THIS FILE IS UNDER RCS - DO NOT MODIFY UNLESS YOU HAVE
#   CHECKED IT OUT USING THE COMMAND CHECKOUT.
#
#

CFLAGS = $(GLOBALFLAGS) -O3 -g -I../include -I/home/benyang/ew/pgplot -flto -Wall
LIBS = -lm -lpng -lz -lX11
FCOMPL = gfortran

L = ../lib

LLIBS = $L/okada91.o $L/matrix.o $L/faultgeom.o

all: hw5_1

hw5_1: hw5_1.o $(LIBS) $(LLIBS)
	$(CC) $(CFLAGS) -o hw5_1 hw5_1.o $(LLIBS) $(LIBS)

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f ./hw5_1
