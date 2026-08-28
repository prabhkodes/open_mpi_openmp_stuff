TARGET = mini_esm
all: $(TARGET)

include ../common.mk

NP  = 8
OMP = 2

# Dependency order matters for Fortran modules. Listed explicitly here so
# you can SEE the chain -- Exercise 21 is about letting CMake infer it, this
# one is about knowing what it inferred.
OBJS = kinds.o mesh_mod.o field_mod.o kernels_mod.o remap_mod.o couple_mod.o mini_esm.o

VPATH = src

$(TARGET): $(OBJS)
	$(MPIFC) $(FFLAGS) $^ -o $@ $(LDLIBS)

%.o: src/%.f90
	$(MPIFC) $(FFLAGS) -c $< -o $@

# Explicit module dependencies (make cannot scan `use` statements).
mesh_mod.o    : kinds.o
field_mod.o   : kinds.o
kernels_mod.o : kinds.o field_mod.o
remap_mod.o   : kinds.o
couple_mod.o  : kinds.o
mini_esm.o    : kinds.o mesh_mod.o field_mod.o kernels_mod.o remap_mod.o couple_mod.o

.PHONY: run spec scaling test
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

spec:
	@cat SPEC.md

## scaling: reuse the Exercise 13 methodology once the model actually runs
scaling: $(TARGET)
	@for p in 2 4 8 16; do \
	  echo "---- $$p ranks ----"; \
	  $(RUNENV) $(MPIRUN) -n $$p ./$(TARGET) 100 64; \
	done

## test: your regression suite (SPEC section 6, built as in Exercise 22)
test: $(TARGET)
	@test -x run_tests.sh && ./run_tests.sh || \
	  echo "run_tests.sh not written yet -- SPEC.md section 6"
