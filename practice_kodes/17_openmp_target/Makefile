TARGET = omp_target
all: $(TARGET)

include ../common.mk

# Host build (macOS / any CPU): target regions execute on the host.
FFLAGS = -O3 -march=native -g -fopenmp -Wall -ffree-line-length-none

$(TARGET): omp_target.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run gpu nvidia sizes
run: $(TARGET)
	$(RUNENV) ./$(TARGET) 192 50

## gpu: GCC with NVPTX offload (Linux cluster with an offload-enabled GCC)
gpu: omp_target.f90
	gfortran -O3 -fopenmp -foffload=nvptx-none -ffree-line-length-none \
	  $< -o omp_target_gpu
	./omp_target_gpu 192 50

## nvidia: NVIDIA HPC SDK -- the toolchain an ESM would actually use
nvidia: omp_target.f90
	nvfortran -O3 -mp=gpu -gpu=cc80 -Minfo=mp $< -o omp_target_nv
	./omp_target_nv 192 50

## sizes: transfer cost dominates more as the array grows
sizes: $(TARGET)
	@for n in 96 128 160 192; do \
	  echo "---- n = $$n ----"; \
	  $(RUNENV) ./$(TARGET) $$n 30 | grep -E 'variant|A:|B:|C:|host'; \
	done

clean-extra:
	rm -f omp_target_gpu omp_target_nv
