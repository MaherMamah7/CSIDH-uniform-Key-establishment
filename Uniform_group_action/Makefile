# Uniform group action: sample a uniform ideal class and act with it on the
# qt-Pegasis starting curve. Compiles the needed qt-Pegasis sources directly,
# so it is self-contained per PNAME (500/1000/1500/2000/4000).
CC = gcc
PNAME ?= 500

GMP_PREFIX ?= $(shell if [ -f /opt/homebrew/include/gmp.h ]; then echo /opt/homebrew; elif [ -f /usr/local/include/gmp.h ]; then echo /usr/local; fi)
ifneq ($(strip $(GMP_PREFIX)),)
INC_GMP = -I$(GMP_PREFIX)/include
LIB_GMP = -L$(GMP_PREFIX)/lib
endif

ROOT = ..
INCLUDES = -I$(ROOT)/src/action/include -I$(ROOT)/src/qt_interface/include -I$(ROOT)/include \
           -I$(ROOT)/src/dim4/include -I$(ROOT)/src/ec/include -I$(ROOT)/src/gf/include \
           -I$(ROOT)/src/gf/test -I$(ROOT)/src/mp/include -I$(ROOT)/src/params \
           -I$(ROOT)/src/common/include -I$(ROOT)/src/quadratic/include \
           -I$(ROOT)/Cayley_uniform_ideal $(INC_GMP)

CFLAGS = -O2 -DNDEBUG -mcpu=native -DP_$(PNAME) $(INCLUDES) -Wno-unused-command-line-argument
LDFLAGS = $(LIB_GMP) -lgmp -lm

SRC = uniform_action.c \
      $(ROOT)/Cayley_uniform_ideal/qform.c \
      $(ROOT)/src/action/qt_pegasis.c \
      $(ROOT)/src/qt_interface/qt_normeq_solve.c \
      $(ROOT)/src/qt_interface/qt_group.c \
      $(ROOT)/src/qt_interface/qt_sos_consts.c \
      $(ROOT)/src/dim4/gluing_isogeny_dim4.c \
      $(ROOT)/src/dim4/isogeny_chain_dim4.c \
      $(ROOT)/src/dim4/isogeny_dim4.c \
      $(ROOT)/src/dim4/splitting_dim4.c \
      $(ROOT)/src/dim4/theta_structures.c \
      $(ROOT)/src/ec/ec.c \
      $(ROOT)/src/ec/basis.c \
      $(ROOT)/src/gf/fp.c \
      $(ROOT)/src/gf/fp/fp_$(PNAME).c \
      $(ROOT)/src/mp/mp.c \
      $(ROOT)/src/common/matrices_mod4.c \
      $(ROOT)/src/common/randombytes_system.c \
      $(ROOT)/src/quadratic/intbig.c \
      $(ROOT)/src/params/fp_params/constants_$(PNAME).c

uniform_action: $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o uniform_action_$(PNAME) $(LDFLAGS)

run: uniform_action
	./uniform_action_$(PNAME)

clean:
	rm -f uniform_action_*
