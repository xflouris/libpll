# Feb 17 2015, converted fortran to C
# Stephen Becker, stephen.becker@colorado.edu

CC = gcc
#CC = gcc-4.9

CFLAGS = -O3 -Wall -Wno-uninitialized 
LDFLAGS = -lpll -lm

INCLUDES = -I/usr/local/include -I./

DRIVER1_77 = driver1.c
DRIVER2_77 = driver2.c
DRIVER3_77 = driver3.c

LBFGSB  = lbfgsb.c linesearch.c subalgorithms.c print.c optimize.c

LINPACK = linpack.c

BLAS 	= miniCBLAS.c

TIMER   = timer.c

SRC = $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

default: all
all :  lbfgsb_pll lbfgsb_77_1 

lbfgsb_77_1 : $(DRIVER1_77C) $(SRC)
	$(CC) $(CFLAGS) $(DRIVER1_77) $(SRC) $(LDFLAGS) -o example_driver1

lbfgsb_pll : $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) pll-lbfgsb.c $(SRC) $(LDFLAGS) -o pll-lbfgsb

lbfgsb_77_3 : $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) $(DRIVER3_77) $(SRC) $(LDFLAGS) -o example_driver3

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
