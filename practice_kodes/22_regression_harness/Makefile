TARGET = model
all: $(TARGET)

include ../common.mk

$(TARGET): model.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run test restartsweep clean-out
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET) run 1000 0 state.txt

## test: the CI suite. Expect test 1 to FAIL before you fix the bug.
test: $(TARGET)
	@chmod +x run_tests.sh
	@./run_tests.sh

## restartsweep: restart at many points -- some expose the bug, some do not
## (TODO 5d). A bug that hides at round numbers is the dangerous kind.
restartsweep: $(TARGET)
	@for r in 100 250 500 617 750 900; do \
	  $(MPIRUN) -n 4 ./$(TARGET) run     1000 0   c.txt >/dev/null 2>&1; \
	  $(MPIRUN) -n 4 ./$(TARGET) restart 1000 $$r r.txt >/dev/null 2>&1; \
	  if cmp -s c.txt r.txt; then res="identical"; else res="DIFFERS"; fi; \
	  printf "  restart at %4s : %s\n" "$$r" "$$res"; \
	done

clean-out:
	rm -f state.txt cont.txt rst.txt c.txt r.txt restart.bin.* rank*.txt
