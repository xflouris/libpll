# Copyright (C) 2015 Tomas Flouri, Diego Darriba
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Profiling
#PROFILING=-g -pg
PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare

AR = ar
CC = gcc
CFLAGS = -D_GNU_SOURCE -fPIC $(PROFILING) -msse3 -O3 -I../ $(WARN)
LINKFLAGS=$(PROFILING)
LIBS=-lpll -lm

PROG=libpll_optimize.so

all: $(PROG)

OBJS=pll_optimize.o \
  opt_algorithms.o \
	lbfgsb/lbfgsb.o lbfgsb/linesearch.o lbfgsb/linpack.o \
	lbfgsb/miniCBLAS.o lbfgsb/subalgorithms.o

$(PROG): $(OBJS)
	$(CC) -shared -Wl,-soname,$@ -o $@ $+

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
