
CC = gcc
BUILD_DIRECTORY = build/
PNAME=500
FLAGS=RELEASE_FLAGS

# Auto-detect a Homebrew/local GMP install (macOS) and add its paths so the
# project builds with a plain `make`. On Linux with a system GMP this stays
# empty and nothing changes.
GMP_PREFIX ?= $(shell if [ -f /opt/homebrew/include/gmp.h ]; then echo /opt/homebrew; elif [ -f /usr/local/include/gmp.h ]; then echo /usr/local; fi)
ifneq ($(strip $(GMP_PREFIX)),)
# Only the include path is needed at compile time. Passing -L on compile-only
# (-c) steps triggers clang's "argument unused during compilation" warning, so
# add -Wno-unused-command-line-argument to keep those steps quiet.
CC := $(CC) -I$(GMP_PREFIX)/include -L$(GMP_PREFIX)/lib -Wno-unused-command-line-argument
endif

#manage the / at end of BUILD_DIRECTORY (which could have been omitted)
TEST = $(patsubst %/,"/",$(BUILD_DIRECTORY))
ifeq ($(TEST),"/")
	BUILDDIR = $(BUILD_DIRECTORY)
else
	BUILDDIR = $(BUILD_DIRECTORY)/
endif

WARNING_FLAGS =  -Wall -Wextra
DISABLE_WARNING_FLAGS =  -Wno-all -Wno-extra
ASAN_FLAGS = -fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1
DEBUG_FLAGS =  -Og -g
RELEASE_FLAGS = -O3 -DNDEBUG -mcpu=native
USE = ${CC} -DP_${PNAME} ${WARNING_FLAGS} ${${FLAGS}}

all: build_all

clean:
	rm -r ${BUILDDIR}*

test: ${BUILDDIR}test_fp_${PNAME} ${BUILDDIR}test_ec_${PNAME} ${BUILDDIR}test_dim4_${PNAME} ${BUILDDIR}test_action_${PNAME} 
	./${BUILDDIR}test_fp_${PNAME}
	./${BUILDDIR}test_ec_${PNAME}
	./${BUILDDIR}test_dim4_${PNAME}
	./${BUILDDIR}test_action_${PNAME}

build_all: ${BUILDDIR}bench_fp_${PNAME} ${BUILDDIR}test_fp_${PNAME} ${BUILDDIR}test_ec_${PNAME} ${BUILDDIR}test_dim4_${PNAME} ${BUILDDIR}test_action_${PNAME} ${BUILDDIR}bench_dim4_${PNAME} 
	@echo "Built all tests and benchmarks\n"

bench: ${BUILDDIR}bench_fp_${PNAME} ${BUILDDIR}bench_dim4_${PNAME} 
	./${BUILDDIR}bench_fp_${PNAME}
	./${BUILDDIR}bench_dim4_${PNAME}

# Full-C end-to-end demo: sample an ideal, solve the norm equation (Step 1) in
# C, then run the 4D action (Steps 2-3), and verify the group-law round trip.
run: ${BUILDDIR}qt_run_${PNAME}
	./${BUILDDIR}qt_run_${PNAME}

# Full-C timing: average per-action time of Step 1 (norm eq) + Steps 2-3 (4D).
bench_full: ${BUILDDIR}qt_bench_${PNAME}
	./${BUILDDIR}qt_bench_${PNAME}

test_arith: ${BUILDDIR}test_fp_${PNAME} 
	./${BUILDDIR}test_fp_${PNAME}

bench_arith: ${BUILDDIR}bench_fp_${PNAME} 
	./${BUILDDIR}bench_fp_${PNAME}

test_theta: ${BUILDDIR}test_fp_${PNAME} ${BUILDDIR}test_dim4_${PNAME} 
	./${BUILDDIR}test_dim4_${PNAME}

test_qt: ${BUILDDIR}bench_qt_${PNAME}
	./${BUILDDIR}bench_qt_${PNAME} -test

isog_test_${PNAME}: ${BUILDDIR}test_dim4_${PNAME}
	./${BUILDDIR}test_dim4_${PNAME}

action_test_${PNAME}: ${BUILDDIR}test_action_${PNAME}
	./${BUILDDIR}test_action_${PNAME}

bench_dim4_${PNAME}: ${BUILDDIR}bench_dim4_${PNAME}
	./${BUILDDIR}bench_dim4_${PNAME}


#Common library
COMMON_INCLUDES=-Isrc/common/include
COMMON_INCLUDEFILES := $(wildcard src/common/include/*.h)
COMMON_CODEFILES=src/common
COMMON_OBJECTFILES=${BUILDDIR}common

COMMON_SOURCES := $(wildcard $(COMMON_CODEFILES)/*.c)
COMMON_OBJECTS := $(patsubst $(COMMON_CODEFILES)/%.c,$(COMMON_OBJECTFILES)/%.o,$(COMMON_SOURCES))

$(COMMON_OBJECTFILES)/%.o: $(COMMON_CODEFILES)/%.c ${COMMON_INCLUDEFILES}
	mkdir -p $(COMMON_OBJECTFILES)
	${USE} -c ${COMMON_INCLUDES}  $< -o $@

${BUILDDIR}libcommon.a: $(COMMON_OBJECTS)
	ar -rcs ${BUILDDIR}libcommon.a  $(COMMON_OBJECTS)


#Params library
PARAMS_INCLUDES = -Isrc/params ${COMMON_INCLUDES} -Isrc/quadratic/include  -Isrc/gf/include -Isrc/mp/include
PARAMS_NEWINCLUDEFILES := $(wildcard src/params/*.h) 
PARAMS_NEWINCLUDEFILES_SND := $(wildcard src/mp/includes/*.h) 
PARAMS_INCLUDEFILES = ${COMMON_INCLUDEFILES} ${PARAMS_NEWINCLUDEFILES} ${PARAMS_NEWINCLUDEFILES_SND}
PARAMS_CODEFILES = src/params/fp_params/constants_${PNAME}.c
PARAMS_OBJECTFILES = ${BUILDDIR}params


${BUILDDIR}params/constants_${PNAME}.o: src/params/fp_params/constants_${PNAME}.c ${PARAMS_INCLUDEFILES}
	mkdir -p $(PARAMS_OBJECTFILES)
	${USE} -c ${PARAMS_INCLUDES}  $< -o $@

${BUILDDIR}libparams_${PNAME}.a: ${BUILDDIR}params/constants_${PNAME}.o
	ar -rcs ${BUILDDIR}libparams_${PNAME}.a  ${BUILDDIR}params/constants_${PNAME}.o

#FP-library
GF_INCLUDES = -Isrc/gf/include -Isrc/mp/include ${PARAMS_INCLUDES}
GF_NEWINCLUDEFILES := $(wildcard src/gf/include/*.h)
GF_INCLUDEFILES = ${PARAMS_INCLUDEFILES} ${MP_INCLUDEFILES} ${GF_NEWINCLUDEFILES}
GF_CODEFILES = src/gf
GF_OBJECTFILES = ${BUILDDIR}gf
GF_SOURCES := $(wildcard $(GF_CODEFILES)/*.c)
GF_OBJECTS := $(patsubst $(GF_CODEFILES)/%.c,$(GF_OBJECTFILES)/%.o,$(GF_SOURCES))

$(GF_OBJECTFILES)/%.o: $(GF_CODEFILES)/%.c ${GF_INCLUDEFILES}
	mkdir -p $(GF_OBJECTFILES)/test
	${USE} -c ${GF_INCLUDES}  $< -o $@

${BUILDDIR}gf/fp_${PNAME}.o: src/gf/fp/fp_${PNAME}.c ${GF_INCLUDEFILES}
	mkdir -p $(GF_OBJECTFILES)
	${USE} ${DISABLE_WARNING_FLAGS} -c ${GF_INCLUDES}  src/gf/fp/fp_${PNAME}.c -o ${BUILDDIR}gf/fp_${PNAME}.o

${BUILDDIR}libgf_${PNAME}.a: $(GF_OBJECTS) ${BUILDDIR}gf/fp_${PNAME}.o
	ar -rcs ${BUILDDIR}libgf_${PNAME}.a  $(GF_OBJECTS) ${BUILDDIR}gf/fp_${PNAME}.o

#FP-tests
GF_TEST_INCLUDES = -Isrc/gf/test -Isrc/mp/include ${GF_INCLUDES}
GF_TEST_NEWINCLUDEFILES := $(wildcard src/gf/test/*.h)
GF_TEST_INCLUDEFILES = ${GF_INCLUDEFILES} ${GF_TEST_NEWINCLUDEFILES}
GF_TEST_CODEFILES = src/gf/test
GF_TEST_OBJECTFILES = ${BUILDDIR}gf/test
GF_TEST_LINK = -lgmp -lm -L./${BUILDDIR} -lgf_${PNAME} -lgf_test_utils_${PNAME} -lparams_${PNAME} -lquadratic -lcommon 
GF_TEST_USED_LIBS = ${BUILDDIR}libgf_${PNAME}.a ${BUILDDIR}libgf_test_utils_${PNAME}.a ${BUILDDIR}libparams_${PNAME}.a ${BUILDDIR}libcommon.a ${BUILDDIR}libquadratic.a 

GF_TEST_SOURCES := $(wildcard $(GF_TEST_CODEFILES)/*.c)
GF_TEST_OBJECTS := $(patsubst $(GF_TEST_CODEFILES)/%.c,$(GF_TEST_OBJECTFILES)/%.o,$(GF_TEST_SOURCES))

$(GF_TEST_OBJECTFILES)/%.o: $(GF_TEST_CODEFILES)/%.c  ${GF_TEST_INCLUDEFILES}
	mkdir -p $(GF_TEST_OBJECTFILES)
	${USE} -c ${GF_TEST_INCLUDES}  $< -o $@

${BUILDDIR}libgf_test_utils_${PNAME}.a: ${BUILDDIR}gf/test/test_utils.o
	ar -rcs ${BUILDDIR}libgf_test_utils_${PNAME}.a  ${BUILDDIR}gf/test/test_utils.o

${BUILDDIR}test_fp_${PNAME}: ${BUILDDIR}gf/test/test_fp.o ${GF_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}test_fp_${PNAME} ${BUILDDIR}gf/test/test_fp.o ${GF_TEST_LINK}
${BUILDDIR}bench_fp_${PNAME}: ${BUILDDIR}gf/test/bench_fp.o ${GF_TEST_USED_LIBS} 
	${USE} -o  ${BUILDDIR}bench_fp_${PNAME} ${BUILDDIR}gf/test/bench_fp.o ${GF_TEST_LINK}

# mp library
MP_INCLUDES=-Isrc/mp/include -Isrc/common/include
MP_NEWINCLUDEFILES := $(wildcard src/mp/include/*.h)
MP_INCLUDEFILES = ${MP_NEWINCLUDEFILES}
MP_CODEFILES=src/mp
MP_OBJECTFILES=${BUILDDIR}mp

MP_SOURCES := $(wildcard $(MP_CODEFILES)/*.c)
MP_OBJECTS := $(patsubst $(MP_CODEFILES)/%.c,$(MP_OBJECTFILES)/%.o,$(MP_SOURCES))

$(MP_OBJECTFILES)/%.o: $(MP_CODEFILES)/%.c ${MP_INCLUDEFILES}
	mkdir -p $(MP_OBJECTFILES)
	${USE} -c ${MP_INCLUDES}  $< -o $@

${BUILDDIR}libmp.a: $(MP_OBJECTS)
	ar -rcs ${BUILDDIR}libmp.a  $(MP_OBJECTS)

# ec library
EC_INCLUDES=-Isrc/ec/include -Isrc/mp/include ${GF_INCLUDES}
EC_NEWINCLUDEFILES := $(wildcard src/ec/include/*.h)
EC_INCLUDEFILES = ${GF_INCLUDEFILES} ${EC_NEWINCLUDEFILES}
EC_CODEFILES=src/ec
EC_OBJECTFILES=${BUILDDIR}ec

EC_SOURCES := $(wildcard $(EC_CODEFILES)/*.c)
EC_OBJECTS := $(patsubst $(EC_CODEFILES)/%.c,$(EC_OBJECTFILES)/%.o,$(EC_SOURCES))

$(EC_OBJECTFILES)/%.o: $(EC_CODEFILES)/%.c ${EC_INCLUDEFILES}
	mkdir -p $(EC_OBJECTFILES)/test
	${USE} -c ${EC_INCLUDES}  $< -o $@

${BUILDDIR}libec.a: $(EC_OBJECTS)
	ar -rcs ${BUILDDIR}libec.a  $(EC_OBJECTS)

# ec test
EC_TEST_INCLUDES = -Isrc/ec/test ${EC_INCLUDES}
EC_TEST_NEWINCLUDEFILES := $(wildcard src/ec/test/*.h)
EC_TEST_NEWINCLUDEFILES_NEXT := $(wildcard src/gf/test/*.h)
EC_TEST_INCLUDEFILES = ${HD_INCLUDEFILES} ${HD_TEST_NEWINCLUDEFILES} ${HD_TEST_NEWINCLUDEFILES_NEXT}
EC_TEST_CODEFILES = src/ec/test
EC_TEST_OBJECTFILES = ${BUILDDIR}ec/test
EC_TEST_LINK = -lm -L./${BUILDDIR} -lec -lgf_test_utils_${PNAME} -lgf_${PNAME} -lparams_${PNAME} -lquadratic -lcommon -lmp -lgmp
EC_TEST_USED_LIBS = ${BUILDDIR}libec.a ${BUILDDIR}libgf_${PNAME}.a ${BUILDDIR}libgf_test_utils_${PNAME}.a ${BUILDDIR}libparams_${PNAME}.a ${BUILDDIR}libcommon.a ${BUILDDIR}libquadratic.a ${BUILDDIR}libmp.a

EC_TEST_SOURCES := $(wildcard $(EC_TEST_CODEFILES)/*.c)
EC_TEST_OBJECTS := $(patsubst $(EC_TEST_CODEFILES)/%.c,$(EC_TEST_OBJECTFILES)/%.o,$(EC_TEST_SOURCES))

$(EC_TEST_OBJECTFILES)/%.o: $(EC_TEST_CODEFILES)/%.c ${EC_TEST_INCLUDEFILES} 
	mkdir -p $(EC_TEST_OBJECTFILES)
	${USE} -c ${EC_TEST_INCLUDES}  $< -o $@

${BUILDDIR}test_ec_${PNAME}: ${BUILDDIR}ec/test/ec_test.o ${BUILDDIR}ec/test/test_extras.o ${EC_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}test_ec_${PNAME} ${BUILDDIR}ec/test/ec_test.o ${BUILDDIR}ec/test/test_extras.o ${EC_TEST_LINK}

test_ec_${PNAME}: ${BUILDDIR}test_ec_${PNAME}
	./${BUILDDIR}test_ec_${PNAME}

# dim4 library
HD_INCLUDES=-Isrc/dim4/include -Isrc/gf/test ${EC_INCLUDES}
HD_NEWINCLUDEFILES := $(wildcard src/dim4/include/*.h)
HD_INCLUDEFILES = ${EC_INCLUDEFILES} ${HD_NEWINCLUDEFILES}
HD_CODEFILES=src/dim4
HD_OBJECTFILES=${BUILDDIR}dim4

HD_SOURCES := $(wildcard $(HD_CODEFILES)/*.c)
HD_OBJECTS := $(patsubst $(HD_CODEFILES)/%.c,$(HD_OBJECTFILES)/%.o,$(HD_SOURCES))

$(HD_OBJECTFILES)/%.o: $(HD_CODEFILES)/%.c ${HD_INCLUDEFILES}
	mkdir -p $(HD_OBJECTFILES)/test
	${USE} -c ${HD_INCLUDES}  $< -o $@

${BUILDDIR}libdim4.a: $(HD_OBJECTS)
	ar -rcs ${BUILDDIR}libdim4.a  $(HD_OBJECTS)

# dim 4 test
HD_TEST_INCLUDES = -Isrc/gf/test ${HD_INCLUDES}
HD_TEST_NEWINCLUDEFILES := $(wildcard src/dim4/test/*.h)
HD_TEST_NEWINCLUDEFILES_NEXT := $(wildcard src/gf/test/*.h)
HD_TEST_INCLUDEFILES = ${HD_INCLUDEFILES} ${HD_TEST_NEWINCLUDEFILES} ${HD_TEST_NEWINCLUDEFILES_NEXT}
HD_TEST_CODEFILES = src/dim4/test
HD_TEST_OBJECTFILES = ${BUILDDIR}dim4/test
HD_TEST_LINK = -lgmp -lm -L./${BUILDDIR} -ldim4 -lec -lgf_test_utils_${PNAME} -lgf_${PNAME} -lparams_${PNAME} -lquadratic -lcommon  -lmp 
HD_TEST_USED_LIBS = ${BUILDDIR}libdim4.a ${BUILDDIR}libec.a ${BUILDDIR}libgf_${PNAME}.a ${BUILDDIR}libgf_test_utils_${PNAME}.a ${BUILDDIR}libparams_${PNAME}.a ${BUILDDIR}libcommon.a ${BUILDDIR}libmp.a ${BUILDDIR}libquadratic.a 
HD_TEST_SOURCES := $(wildcard $(HD_TEST_CODEFILES)/*.c)
HD_TEST_OBJECTS := $(patsubst $(HD_TEST_CODEFILES)/%.c,$(HD_TEST_OBJECTFILES)/%.o,$(HD_TEST_SOURCES))

$(HD_TEST_OBJECTFILES)/%.o: $(HD_TEST_CODEFILES)/%.c ${HD_TEST_INCLUDEFILES}
	mkdir -p $(HD_TEST_OBJECTFILES)
	${USE} -c ${HD_TEST_INCLUDES}  $< -o $@

${BUILDDIR}test_dim4_${PNAME}: ${BUILDDIR}dim4/test/test_dim4.o ${HD_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}test_dim4_${PNAME} ${BUILDDIR}dim4/test/test_dim4.o ${HD_TEST_LINK}




#Quadratic library
QUADRATIC_INCLUDES=-Isrc/quadratic/include ${PARAMS_INCLUDES} 
QUADRATIC_NEWINCLUDEFILES := $(wildcard src/quadratic/include/*.h)
QUADRATIC_INCLUDEFILES = ${PARAMS_INCLUDEFILES} ${QUADRATIC_NEWINCLUDEFILES} ${QUADRATIC_NEWINCLUDEFILES_SND}
QUADRATIC_CODEFILES=src/quadratic
QUADRATIC_OBJECTFILES=${BUILDDIR}quadratic

QUADRATIC_SOURCES := $(wildcard $(QUADRATIC_CODEFILES)/*.c)
QUADRATIC_OBJECTS := $(patsubst $(QUADRATIC_CODEFILES)/%.c,$(QUADRATIC_OBJECTFILES)/%.o,$(QUADRATIC_SOURCES))

$(QUADRATIC_OBJECTFILES)/%.o: $(QUADRATIC_CODEFILES)/%.c ${QUADRATIC_INCLUDEFILES}
	mkdir -p $(QUADRATIC_OBJECTFILES)/test
	${USE} -c ${QUADRATIC_INCLUDES}  $< -o $@

${BUILDDIR}libquadratic.a: $(QUADRATIC_OBJECTS)
	ar -rcs ${BUILDDIR}libquadratic.a  $(QUADRATIC_OBJECTS)


#qt_normeq
QT_INCLUDES= ${QUADRATIC_INCLUDES} -Isrc/qt_interface/include -Iinclude
QT_NEWINCLUDEFILES := $(wildcard src/qt_interface/include/*.h)
QT_NEWINCLUDEFILES_SND = $(wildcard include/*.h)
QT_INCLUDEFILES = ${QUADRATIC_INCLUDEFILES} ${QT_NEWINCLUDEFILES} ${QT_NEWINCLUDEFILES_SND}
QT_CODEFILES=src/qt_interface
QT_OBJECTFILES=${BUILDDIR}qt_interface
QT_SOURCES := $(wildcard $(QT_CODEFILES)/*.c)
QT_OBJECTS := $(patsubst $(QT_CODEFILES)/%.c,$(QT_OBJECTFILES)/%.o,$(QT_SOURCES))

$(QT_OBJECTFILES)/%.o: $(QT_CODEFILES)/%.c ${QT_INCLUDEFILES}
	mkdir -p $(QT_OBJECTFILES)/test
	${USE} -c ${QT_INCLUDES}  $< -o $@

${BUILDDIR}libqt.a: $(QT_OBJECTS)
	ar -rcs ${BUILDDIR}libqt.a  $(QT_OBJECTS)


# action library
ACTION_INCLUDES=-Isrc/action/include -Isrc/params/include ${HD_INCLUDES} ${QT_INCLUDES}
ACTION_NEWINCLUDEFILES := $(wildcard src/action/include/*.h)
ACTION_INCLUDEFILES = ${EC_INCLUDEFILES} ${HD_INCLUDEFILES} ${ACTION_NEWINCLUDEFILES}
ACTION_CODEFILES=src/action
ACTION_OBJECTFILES=${BUILDDIR}action

ACTION_SOURCES := $(wildcard $(ACTION_CODEFILES)/*.c)
ACTION_OBJECTS := $(patsubst $(ACTION_CODEFILES)/%.c,$(ACTION_OBJECTFILES)/%.o,$(ACTION_SOURCES))

$(ACTION_OBJECTFILES)/%.o: $(ACTION_CODEFILES)/%.c ${ACTION_INCLUDEFILES}
	mkdir -p ${ACTION_OBJECTFILES}/test
	${USE} -c ${ACTION_INCLUDES}  $< -o $@

${BUILDDIR}libaction.a: $(ACTION_OBJECTS)
	ar -rcs ${BUILDDIR}libaction.a  $(ACTION_OBJECTS)

# action test
ACTION_TEST_INCLUDES = -Isrc/action/include ${HD_INCLUDES} ${QT_INCLUDES}
ACTION_TEST_NEWINCLUDEFILES := $(wildcard src/action/test/*.h)
ACTION_TEST_NEWINCLUDEFILES_NEXT := $(wildcard src/gf/test/*.h)
ACTION_TEST_INCLUDEFILES = ${ACTION_INCLUDEFILES} ${ACTION_TEST_NEWINCLUDEFILES} ${ACTION_TEST_NEWINCLUDEFILES_NEXT} ${QT_INCLUDEFILES}
ACTION_TEST_CODEFILES = src/action/test
ACTION_TEST_OBJECTFILES = ${BUILDDIR}action/test
ACTION_TEST_LINK = -lgmp -lm -L./${BUILDDIR} -laction -lqt -ldim4 -lec -lgf_test_utils_${PNAME} -lgf_${PNAME} -lparams_${PNAME} -lquadratic -lcommon -lmp 
ACTION_TEST_USED_LIBS = ${BUILDDIR}libaction.a ${BUILDDIR}libdim4.a ${BUILDDIR}libec.a ${BUILDDIR}libgf_${PNAME}.a ${BUILDDIR}libgf_test_utils_${PNAME}.a ${BUILDDIR}libparams_${PNAME}.a ${BUILDDIR}libcommon.a ${BUILDDIR}libquadratic.a ${BUILDDIR}libmp.a ${BUILDDIR}libqt.a

ACTION_TEST_SOURCES := $(wildcard $(ACTION_TEST_CODEFILES)/*.c)
ACTION_TEST_OBJECTS := $(patsubst $(ACTION_TEST_CODEFILES)/%.c,$(ACTION_TEST_OBJECTFILES)/%.o,$(ACTION_TEST_SOURCES))

$(ACTION_TEST_OBJECTFILES)/%.o: $(ACTION_TEST_CODEFILES)/%.c ${ACTION_TEST_INCLUDEFILES}
	mkdir -p ${ACTION_TEST_OBJECTFILES}
	${USE} -c ${ACTION_TEST_INCLUDES}  $< -o $@

${BUILDDIR}test_action_${PNAME}: ${BUILDDIR}action/test/qt_pegasis_test.o ${ACTION_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}test_action_${PNAME} ${BUILDDIR}action/test/qt_pegasis_test.o ${ACTION_TEST_LINK}

${BUILDDIR}bench_dim4_${PNAME}: ${BUILDDIR}action/test/qt_pegasis_dim4_bench.o ${ACTION_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}bench_dim4_${PNAME} ${BUILDDIR}action/test/qt_pegasis_dim4_bench.o ${ACTION_TEST_LINK}

${BUILDDIR}qt_run_${PNAME}: ${BUILDDIR}action/test/qt_normeq_e2e.o ${ACTION_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}qt_run_${PNAME} ${BUILDDIR}action/test/qt_normeq_e2e.o ${ACTION_TEST_LINK}

${BUILDDIR}qt_bench_${PNAME}: ${BUILDDIR}action/test/qt_bench_full.o ${ACTION_TEST_USED_LIBS}
	${USE} -o  ${BUILDDIR}qt_bench_${PNAME} ${BUILDDIR}action/test/qt_bench_full.o ${ACTION_TEST_LINK}

