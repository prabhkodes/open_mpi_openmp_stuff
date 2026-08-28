TARGET = gpu_halo
all: $(TARGET)

include ../common.mk

FFLAGS = -O3 -g -fopenacc -fopenmp -Wall -ffree-line-length-none

$(TARGET): gpu_halo.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run nvidia gpucheck tiles
run: $(TARGET)
	$(MPIRUN) -n $(NP) ./$(TARGET)

## nvidia: the real toolchain on a GPU cluster
nvidia: gpu_halo.f90
	OMPI_FC=nvfortran mpifort -O3 -acc -gpu=cc80 -Minfo=accel $< -o gpu_halo_nv
	OMPI_MCA_opal_cuda_support=1 mpirun -n $(NP) ./gpu_halo_nv

## gpucheck: is this MPI actually GPU-aware? Run this BEFORE you debug
## a segfault in variant B (TODO 5b).
gpucheck:
	@echo "--- OpenMPI CUDA support ---"
	@ompi_info --parsable --all 2>/dev/null | grep -i cuda_support || \
	  echo "  (ompi_info found nothing -- not an OpenMPI build, or no CUDA)"
	@echo "--- UCX devices ---"
	@ucx_info -d 2>/dev/null | grep -i -A1 cuda || echo "  (no ucx_info)"
	@echo "--- visible GPUs ---"
	@nvidia-smi -L 2>/dev/null || echo "  (no nvidia-smi -- CPU-only host)"

## tiles: sweep the interior:halo ratio (TODO 5c)
tiles: $(TARGET)
	@for s in "128 128" "256 256" "512 512" "1024 1024"; do \
	  set -- $$s; \
	  echo "---- tile $$1 x $$2 ----"; \
	  $(MPIRUN) -n $(NP) ./$(TARGET) $$1 $$2 100 | grep -E 'A:|B:|C:'; \
	done

clean-extra:
	rm -f gpu_halo_nv
