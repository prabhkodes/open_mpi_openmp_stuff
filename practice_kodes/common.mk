# ---------------------------------------------------------------------------
# Shared build settings for every exercise in practice_kodes/
#
# Each exercise Makefile does:
#     TARGET = foo
#     all: $(TARGET)
#     include ../common.mk
#     $(TARGET): foo.f90
#         $(MPIFC) $(FFLAGS) $< -o $@
#
# Override anything on the command line:
#     make run NP=8 OMP=4
#     make FFLAGS="-O3 -march=native -fopenmp"
# ---------------------------------------------------------------------------

MPIFC  ?= mpifort
MPICC  ?= mpicc
MPICXX ?= mpicxx

# Serial / pure-OpenMP exercises don't need the MPI wrappers.
# NOTE: plain '=' not '?=' — GNU make predefines FC=f77 and CC=cc as built-in
# variables, so '?=' would silently keep those and try to invoke f77.
# A command-line override (make FC=ifx) still wins over a '=' assignment.
FC = gfortran
CC = gcc-15

# OpenMPI's mpicc wraps Apple clang on macOS, which has no built-in OpenMP.
# Route C/C++ through Homebrew GCC so -fopenmp works the same way everywhere.
ifeq ($(shell uname -s),Darwin)
  export OMPI_CC  ?= gcc-15
  export OMPI_CXX ?= g++-15
endif

# -fcheck=bounds is deliberate: these are learning exercises, and an
# out-of-bounds halo index should abort loudly, not silently corrupt a field.
FFLAGS  ?= -O2 -g -fopenmp -Wall -fcheck=bounds -ffree-line-length-none
CFLAGS  ?= -O2 -g -fopenmp -Wall -std=c11
CXXFLAGS?= -O2 -g -fopenmp -Wall -std=c++17
LDLIBS  ?= -lm

NP     ?= 4
OMP    ?= 2
MPIRUN ?= mpirun --oversubscribe

# Thread pinning matters a lot for hybrid runs on a real cluster, but macOS
# has no thread-affinity API, so libgomp warns on every launch. Only set the
# affinity variables where they mean something.
ifeq ($(shell uname -s),Darwin)
  RUNENV = OMP_NUM_THREADS=$(OMP)
else
  RUNENV = OMP_NUM_THREADS=$(OMP) OMP_PROC_BIND=close OMP_PLACES=cores
endif

.PHONY: clean
clean:
	@rm -f $(TARGET) *.o *.mod *.smod
	@rm -rf *.dSYM
