TARGET = omp_reduce
all: $(TARGET)

include ../common.mk

$(TARGET): omp_reduce.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run repro scale
run: $(TARGET)
	$(RUNENV) ./$(TARGET)

## repro: the ordered sum must print the same hex for every thread count
repro: $(TARGET)
	@for t in 1 2 4 8; do \
	  echo "---- OMP_NUM_THREADS=$$t ----"; \
	  OMP_NUM_THREADS=$$t ./$(TARGET) | grep -E 'bit pattern|ordered'; \
	done

## scale: thread scaling of the whole exercise
scale: $(TARGET)
	@for t in 1 2 4 8; do \
	  echo "---- OMP_NUM_THREADS=$$t ----"; \
	  OMP_NUM_THREADS=$$t ./$(TARGET); \
	done
