#include <mpi.h>
#include <iostream>
#include <vector>

#include "../include/gen_sleep_times.hpp"
#include "../include/supervise.hpp"
#include "../include/work.hpp"


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        std::cerr << "World size must be greater than 1 for " << argv[0] << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        std::vector<int> sleep_times;
        generate_sleep_times(sleep_times, size); // assuming 1 job per process
        supervise_work(sleep_times, size);
    } else {
        do_work();
    }

    MPI_Finalize();
    return 0;
}