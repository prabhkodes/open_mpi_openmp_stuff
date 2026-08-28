TARGET = stencil
all: $(TARGET)

include ../common.mk

CFLAGS = -O3 -march=native -g -fopenmp -Wall -std=c11

$(TARGET): stencil.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

.PHONY: run blocksweep threads sizes
run: $(TARGET)
	$(RUNENV) ./$(TARGET) 256 16 16

## blocksweep: find the blocking that keeps the working set in L2 (TODO 6b)
blocksweep: $(TARGET)
	@for bj in 4 8 16 32 64 128; do \
	  for bk in 4 8 16 32; do \
	    printf "bj=%-4s bk=%-4s  " $$bj $$bk; \
	    $(RUNENV) ./$(TARGET) 256 $$bj $$bk | grep blocked | head -1; \
	  done; \
	done

## threads: does the kernel scale with threads, or is it bandwidth-capped?
threads: $(TARGET)
	@for t in 1 2 4 8; do \
	  echo "---- OMP_NUM_THREADS=$$t ----"; \
	  OMP_NUM_THREADS=$$t ./$(TARGET) 256 16 16 | grep -E 'STREAM|naive|blocked'; \
	done

## sizes: cross from cache-resident to DRAM-bound
sizes: $(TARGET)
	@for n in 64 128 192 256 320; do \
	  echo "---- $${n}^3 ----"; \
	  $(RUNENV) ./$(TARGET) $$n 16 16 | grep -E 'per array|naive|blocked'; \
	done
