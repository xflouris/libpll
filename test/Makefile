#    Copyright (C) 2015 Diego Darriba
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
#    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
#    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

CC = gcc
CFLAGS = -L. -g -O3 -Wall -std=c99
CLIBS = -lpll -lm

CFILES = $(shell find src -name '*.c' ! -name 'common.c')

OBJFILES = $(patsubst src/%.c, obj/%, $(CFILES))

OBJCOMMON = src/common.o

DATADIR = testdata
RESULTDIR = result

DEPFILES = worms16s.fas \
					 small.fas small.tree small.rooted.tree small.rooted.tip.tree \
					 246x4465.fas 246x4465.tree \
					 ribosomal_l5_pf00673.fas ribosomal_l5_pf00673.tree \
					 medium.fas medium.tree \
					 2000.fas 2000.tree \
					 200.fas 200.tree \
					 1000x5.dna.fas 1000.tree

REQFILES = $(patsubst %, $(DATADIR)/%, $(DEPFILES))
ASSETS = https://raw.githubusercontent.com/xflouris/assets/master/libpll

all: $(OBJCOMMON) $(OBJFILES) $(REQFILES)
	@mkdir -p $(RESULTDIR)

$(DATADIR)/%:
	@mkdir -p "$(@D)"
	wget -O $@ $(ASSETS)/$@

obj/%: src/%.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(OBJCOMMON) $(CLIBS) $(LDFLAGS) -DDATADIR=\"$(DATADIR)\"

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $< $(CLIBS) $(LDFLAGS) -DDATADIR=\"$(DATADIR)\"
	
clean:
	rm -rf obj result/* $(OBJCOMMON)
