# Feb 17 2015, converted fortran to C
# Stephen Becker, stephen.becker@colorado.edu

CC = gcc
#CC = gcc-4.9

CFLAGS = -g -O3 -Wall -Wno-uninitialized
LDFLAGS = -lpll_optimize -lpll -lm

INCLUDES = -I/usr/local/include -I./

#LBFGSB  = lbfgsb.c linesearch.c subalgorithms.c optimize.c

#LINPACK = linpack.c

#BLAS 	= miniCBLAS.c

#TIMER   = timer.c

OBJS = src/utils.o src/model-optimize.o
SRC = 
PROG = model-optimize
# $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

default: all
all : $(PROG)

model-optimize : $(OBJS) $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) $(OBJS) $(INCLUDES) $(LDFLAGS) -o $(PROG)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
