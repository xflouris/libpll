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

OS = $(shell uname -s)
CC ?= gcc
PROG=libpll.so

CFLAGS = -DHAVE_SSE -DHAVE_AVX -D_GNU_SOURCE -fPIC -std=c99 -g -msse3 -O3 $(WARN)
AVXFLAGS = -mavx
SSEFLAGS = -msse3

ifeq ($(OS), Darwin)
  SONAMEFLAGS=-Wl,-install_name,$(PROG)
else
  SONAMEFLAGS=-Wl,-soname,$(PROG)
endif

LINKFLAGS=$(PROFILING)
LIBS=-lm

BISON = bison
FLEX = flex


all: $(PROG)

OBJS=fasta.o gamma.o likelihood.o list.o maps.o models.o pll.o output.o \
     parse_utree.o parse_rtree.o lex_utree.o lex_rtree.o utree.o rtree.o \
     derivatives.o partials.o compress.o parse_phylip.o lex_phylip.o \
     utree_moves.o parsimony.o core_derivatives.o core_partials.o \
     core_pmatrix.o core_partials_avx.o core_derivatives_avx.o \
     core_likelihood.o core_likelihood_avx.o core_partials_sse.o \
     core_derivatives_sse.o core_likelihood_sse.o core_pmatrix_avx.o

$(PROG): $(OBJS)
	$(CC) -shared $(SONAMEFLAGS) -o $@ $+

%_avx.o: %_avx.c
	$(CC) $(CFLAGS) $(AVXFLAGS) -c -o $@ $<

%_sse.o: %_sse.c
	$(CC) $(CFLAGS) $(SSEFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

parse_%.c: parse_%.y
	$(BISON) -p pll_$*_ -d -o $@ $<

lex_%.c: lex_%.l
	$(FLEX) -P pll_$*_ -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG) parse_rtree.h parse_utree.h parse_phylip.h parse_rtree.c parse_utree.c parse_phylip.c
