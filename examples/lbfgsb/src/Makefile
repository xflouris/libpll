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

SRC = 
# $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

default: all
all :  lbfgsb_pll 

lbfgsb_pll : $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) pll-lbfgsb.c $(INCLUDES) $(LDFLAGS) -o pll-lbfgsb

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
