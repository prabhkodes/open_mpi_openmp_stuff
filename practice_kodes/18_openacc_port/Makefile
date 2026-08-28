TARGET = acc_port
all: $(TARGET)

include ../common.mk

# -fopenacc on host: gfortran compiles the directives and runs them on the
# CPU. -fopt-info-optimized is the poor-relation of -Minfo=accel, but it is
# what you have locally.
FFLAGS = -O3 -march=native -g -fopenacc -fopenmp -Wall -ffree-line-length-none

$(TARGET): acc_port.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run minfo nvidia
run: $(TARGET)
	$(RUNENV) ./$(TARGET) 192 50

## minfo: what did the compiler actually do? Never tune OpenACC without this.
minfo: acc_port.f90
	$(FC) $(FFLAGS) -fopt-info-optimized -c $< -o /dev/null 2>&1 | head -40

## nvidia: the real toolchain. -Minfo=accel is the output that matters --
## read every line before you look at a single timing.
nvidia: acc_port.f90
	nvfortran -O3 -acc -gpu=cc80 -Minfo=accel $< -o acc_port_nv
	./acc_port_nv 192 50

clean-extra:
	rm -f acc_port_nv
