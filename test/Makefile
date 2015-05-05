CC = gcc
CFLAGS = -g -O3 -Wall 
CLIBS = -lpll -lm

CFILES = $(wildcard **/*.c)
CFILES = $(shell find src/ -type f -name '*.c')
OBJFILES = $(patsubst src/%.c, obj/%, $(CFILES))
DEPS = 

all: $(OBJFILES)

obj/%: src/%.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) -o $@ $< $(CLIBS)

clean:
	rm -rf obj
