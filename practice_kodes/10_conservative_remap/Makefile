TARGET = remap
all: $(TARGET)

include ../common.mk

$(TARGET): remap.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run drift resolutions
run: $(TARGET)
	./$(TARGET)

## drift: long round-trip run -- the climate-drift argument (TODO 5b)
drift: $(TARGET)
	./$(TARGET) 90 143 100000

## resolutions: does conservation hold at every grid ratio?
## It must -- conservation is exact, not asymptotic. If any row fails,
## your overlap sweep has a gap or a double-count.
resolutions: $(TARGET)
	@for r in "10 7" "90 143" "1000 337" "337 1000" "4096 4095"; do \
	  set -- $$r; \
	  echo "---- $$1 -> $$2 ----"; \
	  ./$(TARGET) $$1 $$2 10 | grep -E 'conservation|consistency'; \
	done
