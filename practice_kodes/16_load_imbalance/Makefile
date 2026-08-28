TARGET = imbalance
all: $(TARGET)

include ../common.mk

$(TARGET): imbalance.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run ranks profiled
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET)

## ranks: does the imbalance get worse with more ranks? Predict first.
ranks: $(TARGET)
	@for p in 2 4 8 16; do \
	  echo "---- $$p ranks ----"; \
	  $(MPIRUN) -n $$p ./$(TARGET) 4000 10; \
	done

## profiled: run under the Exercise 15 profiler. Its imbalance analysis
## should independently identify the same slow rank (TODO 5c).
profiled: imbalance.f90 ../15_pmpi_profiler/mpiprof.c
	$(MPICC) $(CFLAGS) -c ../15_pmpi_profiler/mpiprof.c -o mpiprof.o
	$(MPIFC) $(FFLAGS) imbalance.f90 mpiprof.o -o imbalance_profiled
	$(MPIRUN) -n $(NP) ./imbalance_profiled

clean-extra:
	rm -f imbalance_profiled mpiprof.o
