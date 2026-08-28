TARGET = nproma
all: $(TARGET)

include ../common.mk

# -O3 -march=native matters here: the whole exercise is about whether the
# inner loop vectorises, and -O2 will not always do it.
FFLAGS = -O3 -march=native -g -fopenmp -Wall -ffree-line-length-none

$(TARGET): nproma.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run permute vecreport threads
run: $(TARGET)
	$(RUNENV) ./$(TARGET)

## permute: same run with randomised cell numbering (TODO 6c)
permute: $(TARGET)
	$(RUNENV) ./$(TARGET) 500000 90 permute

## vecreport: did the inner loop actually vectorise? (TODO 6b)
vecreport: nproma.f90
	$(FC) $(FFLAGS) -fopt-info-vec -c $< -o /dev/null 2>&1 | grep -i "loop vectorized" || \
	  echo "(no loops vectorized -- that is your answer to TODO 6b)"

## threads: thread scaling at fixed nproma (TODO 6e)
threads: $(TARGET)
	@for t in 1 2 4 8; do \
	  echo "---- OMP_NUM_THREADS=$$t ----"; \
	  OMP_NUM_THREADS=$$t ./$(TARGET) 500000 90 | tail -6; \
	done
