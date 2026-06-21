CXX      = mpicxx
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -march=native

# Homebrew dependency paths (macOS)
LIBOMP_PREFIX   := $(shell brew --prefix libomp)
OPENBLAS_PREFIX := $(shell brew --prefix openblas)

INCLUDES = -I./include \
           -I$(OPENBLAS_PREFIX)/include \
           -I$(LIBOMP_PREFIX)/include

LDFLAGS  = -L$(OPENBLAS_PREFIX)/lib -lopenblas \
           -L$(LIBOMP_PREFIX)/lib   -lomp \
           -Wl,-rpath,$(OPENBLAS_PREFIX)/lib \
           -Wl,-rpath,$(LIBOMP_PREFIX)/lib

OMP_FLAGS = -Xpreprocessor -fopenmp

SRC    = src/main.cpp
TARGET = app_mpi_omp.x

# Default: 4 MPI ranks, 2 OpenMP threads each
NP  ?= 4
OMP ?= 2

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $(INCLUDES) $< $(LDFLAGS) -o $@

run: $(TARGET)
	OMP_NUM_THREADS=$(OMP) mpirun -n $(NP) ./$(TARGET)

clean:
	rm -f $(TARGET) statistics.txt
