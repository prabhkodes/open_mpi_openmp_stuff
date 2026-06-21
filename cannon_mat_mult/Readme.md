# Cannon's Algorithm (MPI + OpenMP + BLAS)

Distributed matrix multiplication on a 2D MPI process grid using Cannon's algorithm. Each rank holds a local tile, does an initial skew, then alternates between local DGEMM (via OpenBLAS) and shifting tiles left/up with `MPI_Sendrecv_replace`. OpenMP parallelises the tile initialisation. Per-phase timings are gathered across all ranks and written to `statistics.txt`.

NP must be a perfect square (4, 9, 16, …). N is set to 10 000 in `main.cpp`.

## Build & Run

Requires `open-mpi`, `openblas`, `libomp` (Homebrew on macOS).

```bash
make
make run              # mpirun -n 4, OMP_NUM_THREADS=2
make run NP=9 OMP=4
make clean
```

## Results

Sample — 4 ranks, 2 OMP threads, N=10 000 (Apple M-series):

| Phase | Avg |
|---|---|
| Init A, B | ~200 ms |
| DGEMM per step | ~2.8 s |
| Shift A, B per step | ~15 ms |
| Total | ~9 s |