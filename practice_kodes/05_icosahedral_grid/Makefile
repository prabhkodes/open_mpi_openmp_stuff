TARGET = icogrid
all: $(TARGET)

include ../common.mk

$(TARGET): icogrid.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run levels
run: $(TARGET)
	./$(TARGET) 3

## levels: refinement sweep -- the table for TODO 6a
levels: $(TARGET)
	@for k in 0 1 2 3 4 5 6; do \
	  echo "==== nrefine = $$k ===="; \
	  ./$(TARGET) $$k | tail -n 14; \
	done
