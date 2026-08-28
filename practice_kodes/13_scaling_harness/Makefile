TARGET = scaling
all: $(TARGET)

include ../common.mk

RANKS ?= 1 2 4 8 16

$(TARGET): scaling.f90
	$(MPIFC) $(FFLAGS) $< -o $@

.PHONY: run strong weak plot slurm
run: $(TARGET)
	$(RUNENV) $(MPIRUN) -n $(NP) ./$(TARGET) strong

## strong: fixed global problem, rising rank count -> scaling_strong.csv
strong: $(TARGET)
	@rm -f scaling_strong.csv
	@for p in $(RANKS); do \
	  $(MPIRUN) -n $$p ./$(TARGET) strong 200 >> scaling_strong.csv; \
	done
	@echo "--- scaling_strong.csv ---"; cat scaling_strong.csv

## weak: fixed work PER RANK -> scaling_weak.csv
weak: $(TARGET)
	@rm -f scaling_weak.csv
	@for p in $(RANKS); do \
	  $(MPIRUN) -n $$p ./$(TARGET) weak 200 >> scaling_weak.csv; \
	done
	@echo "--- scaling_weak.csv ---"; cat scaling_weak.csv

plot: strong weak
	python3 plot_scaling.py scaling_strong.csv scaling_weak.csv -o scaling.png

## slurm: emit a batch script for a real cluster (Levante-shaped).
## Reuse this for every cluster run the rest of the week.
slurm:
	@printf '%s\n' \
	  '#!/bin/bash' \
	  '#SBATCH --job-name=scaling' \
	  '#SBATCH --account=CHANGEME' \
	  '#SBATCH --partition=compute' \
	  '#SBATCH --nodes=8' \
	  '#SBATCH --ntasks-per-node=128' \
	  '#SBATCH --time=00:30:00' \
	  '#SBATCH --output=scaling_%j.out' \
	  '' \
	  'set -euo pipefail' \
	  'module load openmpi' \
	  '' \
	  '# Strong scaling: fixed global problem.' \
	  'for p in 1 2 4 8 16 32 64 128 256 512 1024; do' \
	  '  srun -n $$p ./scaling strong 200 >> scaling_strong.csv' \
	  'done' \
	  '' \
	  '# Weak scaling: fixed work per rank.' \
	  'for p in 1 2 4 8 16 32 64 128 256 512 1024; do' \
	  '  srun -n $$p ./scaling weak 200 >> scaling_weak.csv' \
	  'done' \
	  > job_scaling.sh
	@chmod +x job_scaling.sh
	@echo "wrote job_scaling.sh -- set --account before submitting"

clean-csv:
	rm -f scaling_strong.csv scaling_weak.csv scaling.png job_scaling.sh
