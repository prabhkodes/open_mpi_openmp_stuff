TARGET = uhalo
all: $(TARGET)

include ../common.mk

$(TARGET): uhalo.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run strong shape
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

## strong: fixed mesh, rising rank count -- watch halo/owned grow (TODO 6a)
strong: $(TARGET)
	@for p in 1 2 4 8 16; do \
	  echo "---- $$p ranks ----"; \
	  $(MPIRUN) -n $$p ./$(TARGET) 400 400 200 | head -8; \
	done

## shape: same cell count, different aspect ratio (TODO 6b)
shape: $(TARGET)
	@for s in "400 400" "1600 100" "100 1600"; do \
	  set -- $$s; \
	  echo "---- mesh $$1 x $$2 ----"; \
	  $(MPIRUN) -n $(NP) ./$(TARGET) $$1 $$2 200 | head -8; \
	done
