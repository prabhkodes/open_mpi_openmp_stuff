TARGET = hybrid
all: $(TARGET)

include ../common.mk

$(TARGET): hybrid.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run ratio split
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

## ratio: sweep the compute-to-communication ratio (TODO 6a)
ratio: $(TARGET)
	@for s in 1 2 5 10 50 200; do \
	  echo "---- nsweeps = $$s ----"; \
	  $(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET) 8000000 $$s 20; \
	done

## split: fixed 8 cores, vary the ranks x threads split (TODO 6c)
split: $(TARGET)
	@for cfg in "8 1" "4 2" "2 4" "1 8"; do \
	  set -- $$cfg; \
	  echo "---- $$1 ranks x $$2 threads ----"; \
	  OMP_NUM_THREADS=$$2 $(MPIRUN) -n $$1 ./$(TARGET); \
	done
