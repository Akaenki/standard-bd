SHELL=/bin/bash

# Compiler, I only test this using gcc/5.3
CC=gcc-5

# Compiler optimization level. (O3 is recommended)
OPT_LEVEL=-O3 #-qopt-report

# OpenMP flag, -fopenmp for gcc, -openmp for icc. Leave it blank if use clang
OMPFLAG=-fopenmp

# Header and lib path for openBLAS, change them to ur install path
HPATH=/Users/Linling/Documents/OpenBLAS/include/
LIBPATH=/Users/Linling/Documents/OpenBLAS/lib


################
# OTHER INPUTS #
################
# Number of chains
NUM_CHAINS=1 
# Box size if use PBC
BOX_SIZE=20.0
# Lenard-Jones parameter
LJ_EPSILON=0.31
# Spring Constant
SPR_CONST=200.0
# Delta T
DELTA_T=0.001
# Total time steps
TS_TOTAL=10000000
# Ewald sum parameter
K_MAX=4
################

# Common program arguments.
COMMON_PROG_ARGS= \
								 -std=c99 \
								 $(OMPFLAG) \
								 $(OPT_LEVEL) \
								 #-lrt \

# Input arguments.
INPUT_ARGS = \
		-DNP=${NUM_CHAINS} \
		-DL=${BOX_SIZE} \
		-DEPSILON=${LJ_EPSILON} \
		-DKAPPA=${SPR_CONST} \
		-DDT=${DELTA_T} \
		-DTMAX=${TS_TOTAL} \
		-DKMAX=${K_MAX} \
								
# Program arguments.
OMP_CC = $(CC) $(COMMON_PROG_ARGS) $(INPUT_ARGS)

SRC=*.c
TARGETS=sc_omp

all: $(TARGETS)

sc_omp: $(SRC)
	$(OMP_CC) $(SRC) -lopenblas -I$(HPATH) -L$(LIBPATH) -o $@
	if [ ! -d ./output ];then \
	mkdir ./output;	\
	fi

clean:
	rm -f *.o *.optrpt $(TARGETS)
