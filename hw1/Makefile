#
#   THIS FILE IS UNDER RCS - DO NOT MODIFY UNLESS YOU HAVE
#   CHECKED IT OUT USING THE COMMAND CHECKOUT.
#
#

CFLAGS = $(GLOBALFLAGS) -O3 -g -I../include -I/home/benyang/ew/pgplot -flto
LIBS = -lm -lpng -lz -lX11
FCOMPL = gfortran

L = ../lib

LLIBS = $L/doppler.o $L/orbit.o

CPGPLOT_LIB = /home/benyang/ew/pgplot/libcpgplot.a /home/benyang/ew/pgplot/libpgplot.a

hw1_1: $(LIBS) $(LLIBS)
	$(CC) $(CFLAGS) -o ./hw1_1 hw1_1.c $(LIBS) $(LLIBS)

hw1_2: hw1_2.o $(LIBS) $(LLIBS)
	$(FCOMPL) -o hw1_2 hw1_2.o $(LLIBS) $(CPGPLOT_LIB) $(LIBS)

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f ./hw1_1* ./hw1_2
