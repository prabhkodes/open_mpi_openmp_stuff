# MPI Matrix with Operator Overloading

Templated `CMatrix<T, N>` class that works in an MPI context. Each rank fills a 100×100 tile and writes it to a file — even ranks write matrix A, odd ranks write matrix B. Uses `operator<<` overloading so the matrix prints directly to any stream. Timings are gathered from all ranks via `MPI_Gatherv` and written to `timings.txt`.

## Build & Run

```bash
make
make run        # mpirun -n 6 by default
make run NP=4   # must be >= 2
make clean
```