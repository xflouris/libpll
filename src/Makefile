# Copyright (C) 2015 Tomas Flouri
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

ifndef CC
CC = gcc
endif
CFLAGS = -DHAVE_AVX -D_GNU_SOURCE -fPIC -std=c99 -g -msse3 -O3 $(WARN)
AVXFLAGS = -mavx
LINKFLAGS=$(PROFILING)
LIBS=-lm

BISON = bison
FLEX = flex

PROG=libpll.so

all: $(PROG)

OBJS=fasta.o gamma.o likelihood.o list.o maps.o models.o pll.o output.o \
     parse_utree.o parse_rtree.o lex_utree.o lex_rtree.o utree.o rtree.o \
     core_likelihood.o compress.o parse_phylip.o lex_phylip.o \
     core_likelihood_avx.o core_pmatrix.o

$(PROG): $(OBJS)
	$(CC) -shared -Wl,-soname,$@ -o $@ $+

%_avx.o: %_avx.c
	$(CC) $(CFLAGS) $(AVXFLAGS) -c -o $@ $<
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

parse_%.c: parse_%.y
	$(BISON) -p pll_$*_ -d -o $@ $<

lex_%.c: lex_%.l
	$(FLEX) -P pll_$*_ -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG) parse_rtree.h parse_utree.h parse_phylip.h parse_rtree.c parse_utree.c parse_phylip.c
