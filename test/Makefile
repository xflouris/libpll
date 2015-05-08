CC = gcc
CFLAGS = -g -O3 -Wall 
CLIBS = -lpll -lm

CFILES = $(wildcard **/*.c)
CFILES = $(shell find src/ -type f -name '*.c')
OBJFILES = $(patsubst src/%.c, obj/%, $(CFILES))

DATADIR = testdata
RESULTDIR = result

DEPFILES = worms16s.fas small.fas
REQFILES = $(patsubst %, $(DATADIR)/%, $(DEPFILES))
ASSETS = https://raw.githubusercontent.com/xflouris/assets/master/libpll

all: $(OBJFILES) $(REQFILES)
	@mkdir -p $(RESULTDIR)

$(DATADIR)/%:
	@mkdir -p "$(@D)"
	wget -O $@ $(ASSETS)/$@

obj/%: src/%.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) -o $@ $< $(CLIBS) -DDATADIR=\"$(DATADIR)\"


clean:
	rm -rf obj
