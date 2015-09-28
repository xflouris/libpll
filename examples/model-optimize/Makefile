# Feb 17 2015, converted fortran to C
# Stephen Becker, stephen.becker@colorado.edu

CC = gcc
#CC = gcc-4.9

CFLAGS = -g -pg -O3 -Wall -Wno-uninitialized
LDFLAGS = -lpll_optimize -lpll -lm

INCLUDES = -I/usr/local/include -I./

#LBFGSB  = lbfgsb.c linesearch.c subalgorithms.c optimize.c

#LINPACK = linpack.c

#BLAS 	= miniCBLAS.c

#TIMER   = timer.c

OBJ = src/utils.o src/model-optimize.o
SRC = 
# $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

default: all
all : model-optimize

model-optimize : $(OBJ) $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) $(OBJ) $(INCLUDES) $(LDFLAGS) -o model-optimize

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
