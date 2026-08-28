TARGET = multirate
all: $(TARGET)

include ../common.mk

$(TARGET): multirate.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run restarts clean-rst
run: $(TARGET)
	./$(TARGET)

## restarts: restart at many different points -- a correct implementation
## reproduces at EVERY one, including points that are not multiples of the
## coupling window. Those are the ones that expose a missing accumulator.
restarts: $(TARGET)
	@for r in 100 250 400 500 617 900; do \
	  echo "---- restart at step $$r ----"; \
	  ./$(TARGET) 1000 $$r | grep -E 'restart is|continuous|restarted'; \
	done

clean-rst:
	rm -f restart.bin
