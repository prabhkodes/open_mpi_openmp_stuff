TARGET = mpmd
all: $(TARGET)

include ../common.mk

NP = 8

$(TARGET): mpmd.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run split
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

## split: sweep the atmosphere:ocean rank ratio (TODO 5 / TODO 6a)
## Predict the optimum from the 3:1 cost ratio BEFORE you run this.
split: $(TARGET)
	@for a in 1 2 3 4 5 6 7; do \
	  echo "---- atm=$$a  oce=$$(( $(NP) - $$a )) ----"; \
	  $(MPIRUN) -n $(NP) ./$(TARGET) $$a 20 5 | grep -E 'ranks|wall|wasted'; \
	done
