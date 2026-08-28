# practice_kodes — DKRZ RSE interview prep

Targeted exercise set for [Research Software Engineer, HPC Code Performance Optimisation](https://dkrz.softgarden.io/job/66326173?l=en) at DKRZ.

**Start here: [PLAN.md](PLAN.md)** — the day-by-day schedule, gap analysis, and exercise index.

**On your phone: [MOBILE.md](MOBILE.md)** — self-contained revision. 29 collapsible flashcards, numbers to quote, and the "what would you optimise?" procedure. No code to run.

## Layout

```
practice_kodes/
├── PLAN.md            # the 6-day plan — read this first
├── common.mk          # shared compiler/run settings, included by every Makefile
├── ring_allreduce.c   # Exercise 00 (pre-existing warm-up)
├── 01_..24_/          # one directory per exercise
└── notes/             # your write-ups — one per day, then INTERVIEW.md
```

Every exercise source file opens with a header comment containing:

- **GOAL** — what you are building
- **WHY (DKRZ)** — why this specific skill appears in the job description
- **TASKS** — numbered `TODO`s in the code
- **ACCEPTANCE** — the number or property you must produce before moving on
- **HINTS** — enough to unstick you without giving away the answer

## Building and running

```bash
cd 01_fortran_mpi_halo_1d && make run NP=4 OMP=2
```

Common overrides:

| Variable | Default | Meaning |
|---|---|---|
| `NP` | 4 | MPI ranks |
| `OMP` | 2 | OpenMP threads per rank |
| `FFLAGS` / `CFLAGS` | `-O2 -g -fopenmp -Wall` | compiler flags |

Skeletons compile as shipped (with the `TODO`s unimplemented) so a build error always means *your* change, never the starting point.

## Toolchain

Verified on macOS with Homebrew. `common.mk` routes C through `gcc-15` on Darwin because OpenMPI's `mpicc` wraps Apple clang, which has no built-in OpenMP.

```bash
brew install open-mpi gcc libomp openblas hdf5-mpi netcdf netcdf-fortran cmake gnuplot
```

Exercises 17–19 need a GPU and an NVIDIA/AMD compiler. Each ships a CPU-host fallback so the code can be written and debugged on macOS, then run for real on a cluster.

## Rule

**Every exercise ends with a number.** A runtime, a speedup, a parallel efficiency, a bandwidth, a conservation error. Record it in `notes/dayN.md`. DKRZ is hiring a performance engineer — the interview goes much better when you can quote your own measurements.
