# MPI + OpenMP Projects

![C++](https://img.shields.io/badge/C++-00599C?style=flat-square&logo=cplusplus&logoColor=white)
![OpenMPI](https://img.shields.io/badge/OpenMPI-364d6e?style=flat-square&logoColor=white)
![OpenMP](https://img.shields.io/badge/OpenMP-006DB8?style=flat-square&logoColor=white)
![OpenBLAS](https://img.shields.io/badge/OpenBLAS-0096D6?style=flat-square&logoColor=white)

Distributed and shared-memory parallel programs using OpenMPI and OpenMP. Covers communication patterns, dynamic scheduling, and hybrid parallelism for HPC workloads.

## Projects

### cannon_mat_mult
Distributed matrix multiplication using Cannon's algorithm on a 2D MPI process grid. Each rank holds a local tile and performs DGEMM via OpenBLAS, with OpenMP parallelising tile initialisation. Tiles are shifted with `MPI_Sendrecv_replace`. Per-phase timings (init, DGEMM, shift) are written to `statistics.txt`. NP must be a perfect square.

```bash
make run              # mpirun -n 4, OMP_NUM_THREADS=2
make run NP=9 OMP=4
```

### job_scheduler
Supervisor/worker dynamic load balancing using MPI point-to-point. Rank 0 generates jobs and dispatches them to the first available worker (`MPI_ANY_SOURCE`). Workers loop until they receive a stop tag. Fast workers automatically pick up more work — no idle time.

```bash
make run        # 1 supervisor + 3 workers
make run NP=8
```

### matrix_overload_vector_mpi
Templated `CMatrix<T, N>` class that works in an MPI context. Even ranks write matrix A, odd ranks write matrix B to separate files. Implements `operator<<` for direct stream output. Timings gathered across all ranks via `MPI_Gatherv`.

```bash
make run        # mpirun -n 6
make run NP=4
```

### mpi_timer
Benchmarks common MPI communication primitives — `MPI_Send`/`Recv`, `MPI_Scatter`, and `MPI_Put` (one-sided) — across `int` and `double` types, comparing bulk vs element-by-element transfers. Timings written to `timings.txt`.

```bash
make run        # mpirun -n 2
make run NP=4
```

## Dependencies

- OpenMPI
- OpenBLAS (cannon_mat_mult)
- libomp

Since I was using mac, and if you are too:
On macOS: `brew install open-mpi openblas libomp`
