# MPI Job Scheduler

Supervisor/worker pattern using MPI point-to-point messaging. Rank 0 generates a list of jobs and dispatches them one at a time to whatever worker finishes first (`MPI_ANY_SOURCE`). Workers loop until they receive a stop tag. Fast workers automatically pick up more work — no idle time.

## Build & Run

```bash
make
make run        # mpirun -n 4 (1 supervisor + 3 workers)
make run NP=8
make clean
```
