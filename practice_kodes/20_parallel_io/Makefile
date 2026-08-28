TARGET = pario
all: $(TARGET)

include ../common.mk

$(TARGET): pario.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run ranks hints hdf5 clean-out
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET)

## ranks: which methods scale? (TODO 6a) -- watch variant D fall over
ranks: $(TARGET)
	@for p in 1 2 4 8 16; do \
	  echo "---- $$p ranks ----"; \
	  $(MPIRUN) -n $$p ./$(TARGET) 2048 2048 | grep -E 'A:|B:|C:|D:'; \
	done

## hints: sweep the collective-buffering aggregator count (TODO 6c)
## Set the hints in write_mpiio_collective, then re-run this.
hints: $(TARGET)
	@for cb in 1 2 4 8; do \
	  echo "---- cb_nodes = $$cb ----"; \
	  ROMIO_HINTS_CB_NODES=$$cb $(MPIRUN) -n $(NP) ./$(TARGET) | grep 'C:'; \
	done

## hdf5: confirm parallel HDF5 is available for TODO 6e
hdf5:
	@h5pcc -show 2>/dev/null || echo "(no h5pcc -- brew install hdf5-mpi)"
	@nf-config --flibs 2>/dev/null || echo "(no netcdf-fortran)"

clean-out:
	rm -f out_rank_*.bin out_indep.bin out_coll.bin out_gather.bin
