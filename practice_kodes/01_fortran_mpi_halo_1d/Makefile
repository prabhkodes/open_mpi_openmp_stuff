TARGET = halo1d
all: $(TARGET)

include ../common.mk

$(TARGET): halo1d.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run sweep
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

## sweep: time all three variants across halo widths 1 -> 4096
sweep: $(TARGET)
	@for h in 1 4 16 64 256 1024 4096; do \
	  echo "---- halo = $$h doubles ($$((h*8)) bytes) ----"; \
	  $(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET) 1000000 $$h 200; \
	done
