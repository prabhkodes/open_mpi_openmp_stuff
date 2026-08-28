TARGET = interop
all: $(TARGET)

include ../common.mk

# Mixed-language link: compile each side with its own wrapper, then link
# with the FORTRAN driver (mpifort) so the Fortran runtime gets pulled in.
# Linking with mpicc instead would fail on missing _gfortran_* symbols.
$(TARGET): interop.o ckernels.o
	$(MPIFC) $(FFLAGS) $^ -o $@ $(LDLIBS)

interop.o: interop.f90
	$(MPIFC) $(FFLAGS) -c $< -o $@

ckernels.o: ckernels.c
	$(MPICC) $(CFLAGS) -c $< -o $@

.PHONY: run syms
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET)

## syms: inspect the symbol names on both sides of the boundary.
## Without bind(C), gfortran would emit `f_scale_` (trailing underscore) and
## the C linker would never find it. This is what bind(C) is really for.
syms: interop.o ckernels.o
	@echo "--- symbols the Fortran object DEFINES or WANTS ---"
	@nm interop.o | grep -iE 'c_axpy|c_column|c_sum_over|c_register|c_grid|f_scale' || true
	@echo "--- symbols the C object DEFINES or WANTS ---"
	@nm ckernels.o | grep -iE 'c_axpy|c_column|c_sum_over|c_register|c_grid|f_scale' || true
