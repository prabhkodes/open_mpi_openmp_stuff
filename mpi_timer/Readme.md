# MPI Communication Timer

Benchmarks common MPI communication patterns — `MPI_Send`/`Recv`, `MPI_Scatter`, and `MPI_Put` (one-sided) — across `int` and `double`, and compares bulk vs element-by-element transfers. Results are gathered from all ranks and written to `timings.txt`.

## Build & Run

```bash
make
make run        # mpirun -n 2 by default
make run NP=4
make clean
```