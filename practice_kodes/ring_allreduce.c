#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * EXERCISE: Implement a ring allreduce without MPI_Allreduce
 * Each rank contributes one value. Final result: sum of all values on all ranks.
 * 
 * Pattern: Rank i sends to (i+1)%nprocs, receives from (i-1+nprocs)%nprocs
 * Phase 1: Reduce phase (move data around ring)
 * Phase 2: Broadcast phase (send results back around)
 */

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    // Each rank starts with its rank number as the value
    double my_value = (double)(rank);
    double sum = my_value;  // Local sum (will accumulate)
    double received_value;
    
    int next_rank = (rank + 1) % nprocs;
    int prev_rank = (rank - 1 + nprocs) % nprocs;
    
    printf("Rank %d: starting value = %.0f\n", rank, my_value);
    
    // TODO: Implement ring reduce (nprocs-1 steps)
    // Each step: send 'sum' to next_rank, receive from prev_rank
    // What MPI calls do you need? (Send/Recv? Isend/Irecv?)
    // Hint: Blocking Send/Recv is simpler first, but watch for deadlock!
    
    for (int step = 0; step < nprocs - 1; step++) {
        // TODO: Send 'sum' to next_rank
        // TODO: Receive into 'received_value' from prev_rank
        // TODO: Add received_value to sum
    }
    
    // Now broadcast: send result back around (another nprocs-1 steps)
    // But now all ranks have the same 'sum', so just send the final result
    sum = 0;  // Reset for second phase
    for (int i = 0; i < nprocs; i++) sum += i;  // Expected final sum
    
    // TODO: Second phase—broadcast the sum to all ranks
    
    printf("Rank %d: final sum = %.0f\n", rank, sum);
    
    MPI_Finalize();
    return 0;
}