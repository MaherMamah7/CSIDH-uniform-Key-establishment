# Near-uniform ideal-class sampler (random walk on the class-group Cayley graph)
CC = gcc
GMP_PREFIX ?= $(shell if [ -f /opt/homebrew/include/gmp.h ]; then echo /opt/homebrew; elif [ -f /usr/local/include/gmp.h ]; then echo /usr/local; fi)
ifneq ($(strip $(GMP_PREFIX)),)
INC = -I$(GMP_PREFIX)/include
LIB = -L$(GMP_PREFIX)/lib
endif
CFLAGS = -O2 -Wall -Wextra $(INC)
LDFLAGS = $(LIB) -lgmp -lm

all: cayley test_qform

cayley: cayley.c qform.c qform.h
	$(CC) $(CFLAGS) cayley.c qform.c -o cayley $(LDFLAGS)

test_qform: test_qform.c qform.c qform.h
	$(CC) $(CFLAGS) test_qform.c qform.c -o test_qform $(LDFLAGS)

test: test_qform
	./test_qform

# Run the sampler at a qt-Pegasis level, e.g.  make run PNAME=500 SEED=1
PNAME ?= 500
SEED ?= 1
run: cayley
	./cayley $(PNAME) -s $(SEED)

clean:
	rm -f cayley test_qform
