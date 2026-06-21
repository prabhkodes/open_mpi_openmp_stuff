CXX = mpic++
CXXFLAGS = -std=c++17 -O2 -g -Wall -Wextra
INCLUDES = -Iinclude
SRC_DIR = src
TARGET = $(SRC_DIR)/app
SRC = $(SRC_DIR)/main.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<

clean:
	rm -f $(TARGET)

run: $(TARGET)
	mpirun -np 2 $(TARGET)
