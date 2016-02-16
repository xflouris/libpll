# Copyright (C) 2015 Diego Darriba 
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
# Contact: Tomas Flouri <Diego.Darriba@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Profiling
#PROFILING=-g -pg
PROFILING=-g

CC = gcc
#CC = gcc-4.9

CFLAGS = -g -O3 -Wall -Wsign-compare
LDFLAGS = -lpll_optimize -lpll -lm

OBJS = tbr.o
PROG = tbr

default: all
all : $(PROG)

tbr: $(OBJS) 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) $(PROG)
