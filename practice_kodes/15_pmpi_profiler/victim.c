/* ==========================================================================
 * victim.c — a test application for the Exercise 15 profiler
 * ==========================================================================
 *
 * DO NOT MODIFY THIS FILE. That is the whole point: your profiler must work
 * on code it has never seen and cannot change.
 *
 * The communication pattern is deliberately simple and fully known, so you
 * can check your profiler against ground truth. It prints exactly what it
 * did just before calling MPI_Finalize -- your profile must agree.
 *
 * Built-in features for you to detect:
 *   - two distinct message sizes (8 B halo, 800 KB bulk) -> two histogram peaks
 *   - rank 1 does 3x the compute of everyone else        -> load imbalance,
 *     and rank 1 should show the LEAST time in collectives
 * ========================================================================== */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NITER    200
#define BULK_N   100000        /* doubles -> 800 KB */

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int left  = (rank - 1 + nprocs) % nprocs;
    int right = (rank + 1) % nprocs;

    double  halo_send = rank, halo_recv = 0.0;
    double *bulk = malloc(BULK_N * sizeof(double));
    for (int i = 0; i < BULK_N; i++) bulk[i] = (double)(i % 13);

    double local = 1.0, global = 0.0;

    for (int it = 0; it < NITER; it++) {

        /* Rank 1 is deliberately slow. A correct imbalance analysis must
         * finger it -- and it will show the LEAST time inside collectives,
         * because it is always the last to arrive. */
        int work = (rank == 1) ? 3 : 1;
        double acc = 0.0;
        for (int w = 0; w < work * 200000; w++) acc += sqrt((double)w);
        if (acc < 0.0) printf(" ");

        /* --- small messages: 8 bytes each, latency-bound --------------- */
        MPI_Sendrecv(&halo_send, 1, MPI_DOUBLE, right, 10,
                     &halo_recv, 1, MPI_DOUBLE, left,  10,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Request req[2];
        MPI_Irecv(&halo_recv, 1, MPI_DOUBLE, left,  20, MPI_COMM_WORLD, &req[0]);
        MPI_Isend(&halo_send, 1, MPI_DOUBLE, right, 20, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);

        /* --- one big message every 10 iterations: bandwidth-bound ------ */
        if (it % 10 == 0) {
            MPI_Bcast(bulk, BULK_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        /* --- collectives: where the waiting shows up ------------------- */
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        long n_bcast = (NITER + 9) / 10;
        printf("\n=== victim.c: ground truth (per rank, %d ranks) ===\n", nprocs);
        printf("  MPI_Sendrecv   %6d calls   %10ld bytes\n",
               NITER, (long)NITER * 8 * 2);
        printf("  MPI_Isend      %6d calls   %10ld bytes\n",
               NITER, (long)NITER * 8);
        printf("  MPI_Irecv      %6d calls   %10ld bytes\n",
               NITER, (long)NITER * 8);
        printf("  MPI_Waitall    %6d calls\n", NITER);
        printf("  MPI_Bcast      %6ld calls   %10ld bytes\n",
               n_bcast, n_bcast * BULK_N * 8);
        printf("  MPI_Allreduce  %6d calls   %10ld bytes\n",
               NITER, (long)NITER * 8);
        printf("  MPI_Barrier    %6d calls\n", NITER);
        printf("  expected slow rank: 1 (3x the compute)\n");
        printf("  Multiply by %d ranks to compare with the profile totals.\n",
               nprocs);
        printf("===================================================\n");
    }

    free(bulk);
    MPI_Finalize();     /* <- your wrapper prints the profile from here */
    return 0;
}
