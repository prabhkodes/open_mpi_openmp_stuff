TARGET = victim_profiled
all: $(TARGET) victim_bare

include ../common.mk

# Link-time interposition: because mpiprof.o defines MPI_Send etc, the
# linker resolves the application's calls to OURS instead of the MPI
# library's. The application source is untouched and unaware.
#
# This is the portable approach. LD_PRELOAD works on Linux, but macOS uses
# DYLD_INSERT_LIBRARIES and restricts it under SIP, so linking directly is
# the method that works everywhere -- including on a cluster where you may
# not be allowed to preload into an MPI launcher.

$(TARGET): victim.o mpiprof.o
	$(MPICC) $(CFLAGS) $^ -o $@ $(LDLIBS)

## victim_bare: the same app WITHOUT the profiler, for overhead comparison
victim_bare: victim.o
	$(MPICC) $(CFLAGS) $^ -o $@ $(LDLIBS)

victim.o: victim.c
	$(MPICC) $(CFLAGS) -c $< -o $@

mpiprof.o: mpiprof.c
	$(MPICC) $(CFLAGS) -c $< -o $@

.PHONY: run overhead syms ex06 ex09
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET)

## overhead: what does the profiler itself cost? (TODO 7e)
overhead: $(TARGET) victim_bare
	@echo "--- without profiler ---"
	@/usr/bin/time -p $(MPIRUN) -n $(NP) ./victim_bare 2>&1 | tail -4
	@echo "--- with profiler ---"
	@/usr/bin/time -p $(MPIRUN) -n $(NP) ./$(TARGET) 2>&1 | tail -4

## syms: prove the interposition actually happened.
## mpiprof.o must DEFINE (T) the MPI_* symbols and REFERENCE (U) the PMPI_*.
syms: mpiprof.o
	@nm mpiprof.o | grep -E ' [TU] _?P?MPI_(Send|Recv|Barrier|Allreduce|Finalize)' | sort

## ex06: profile the unstructured halo exchange with no source changes
ex06: mpiprof.o
	$(MPIFC) $(FFLAGS) ../06_unstructured_halo/uhalo.f90 mpiprof.o -o uhalo_profiled
	$(MPIRUN) -n $(NP) ./uhalo_profiled

## ex09: profile the coupled components -- does it find the slow component?
ex09: mpiprof.o
	$(MPIFC) $(FFLAGS) ../09_mpmd_components/mpmd.f90 mpiprof.o -o mpmd_profiled
	$(MPIRUN) -n 8 ./mpmd_profiled

clean-extra:
	rm -f victim_bare uhalo_profiled mpmd_profiled
