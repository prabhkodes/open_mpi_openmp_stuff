CXX      = mpicxx
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -g
INCLUDES = -I./include

SRC    = src/main.cpp
TARGET = src/app

# Number of MPI ranks (1 supervisor + N-1 workers; must be >= 2)
NP ?= 4

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@

run: $(TARGET)
	mpirun -n $(NP) ./$(TARGET)

clean:
	rm -f $(TARGET)
