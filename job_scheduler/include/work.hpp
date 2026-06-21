#pragma once
#include <mpi.h>
#include <thread>
#include <chrono>

inline void do_work() {
    const int WORKTAG = 0;
    const int RESULTTAG = 1;
    const int STOPTAG = 2;

    const int root_process{0};

    while (true) {
        int recv_buff{0}; 
        MPI_Status status{};

        MPI_Recv(&recv_buff, 1, MPI_INT, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == STOPTAG) break;

        if (status.MPI_TAG == WORKTAG) {
        // sleep if work lol
            std::this_thread::sleep_for(std::chrono::seconds(recv_buff));}
        int result = recv_buff * 2; // just for visibility

        MPI_Send(&result, 1, MPI_INT, root_process, RESULTTAG, MPI_COMM_WORLD);
    }
}
