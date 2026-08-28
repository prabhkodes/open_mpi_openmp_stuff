TARGET = icesheet
all: $(TARGET)

include ../common.mk

NP = 8

$(TARGET): icesheet.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run compare split profiled
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET) A

## compare: the headline result. A vs B on the same rank count.
## Check the total forcing matches BEFORE you believe the speedup.
compare: $(TARGET)
	@echo "############ STRATEGY A: synchronous ############"
	@$(MPIRUN) -n $(NP) ./$(TARGET) A 400
	@echo "############ STRATEGY B: concurrent, lagged ############"
	@$(MPIRUN) -n $(NP) ./$(TARGET) B 400

## split: sweep the ice-sheet rank count (strategy C, TODO 6c).
## Two failure modes at the two ends -- find them.
split: $(TARGET)
	@for ice in 1 2 3 4; do \
	  echo "---- ice ranks = $$ice ----"; \
	  $(MPIRUN) -n $(NP) ./$(TARGET) C 400 $$ice | \
	    grep -E 'wall time|steps per second|late'; \
	done

## profiled: run under the Exercise 15 profiler -- does it see the stall?
profiled: icesheet.f90 ../15_pmpi_profiler/mpiprof.c
	$(MPICC) $(CFLAGS) -c ../15_pmpi_profiler/mpiprof.c -o mpiprof.o
	$(MPIFC) $(FFLAGS) icesheet.f90 mpiprof.o -o icesheet_profiled
	@echo "############ A ############"
	@$(MPIRUN) -n $(NP) ./icesheet_profiled A 400 | tail -25
	@echo "############ B ############"
	@$(MPIRUN) -n $(NP) ./icesheet_profiled B 400 | tail -25

clean-extra:
	rm -f icesheet_profiled mpiprof.o
