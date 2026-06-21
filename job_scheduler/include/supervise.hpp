#pragma once
#include <mpi.h>
#include <vector>
#include <iostream>

void supervise_work(const std::vector<int>& sleep_times, int world_size) {
    const int WORKTAG = 0;
    const int RESULTTAG = 1;
    const int STOPTAG = 2;

    const long int st_size = static_cast<long int>(sleep_times.size());
    int next_st_idx{0};
    int sent_count{0};

    for (int r = 1; r < world_size && next_st_idx < st_size; ++r) {

        std::cout << "[root] send WORK to w=" << r
                  << " val=" << sleep_times[next_st_idx] << " (task#" << next_st_idx << ")\n" << std::flush;

        MPI_Send(&sleep_times[next_st_idx], 1, MPI_INT, r, WORKTAG, MPI_COMM_WORLD);
        ++sent_count;
        ++next_st_idx;
    }

    while (sent_count > 0) { // keep recv 
        int result{0};
        MPI_Status status{};

        MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, RESULTTAG, MPI_COMM_WORLD, &status);

        int w_rank = status.MPI_SOURCE;

        std::cout << "[root] got RESULT from w=" << w_rank << " result=" << result << "\n" << std::flush;

        if (next_st_idx < st_size) {
            std::cout << "[root] send WORK to w=" << w_rank
                      << " val=" << sleep_times[next_st_idx] << " (task#" << next_st_idx << ")\n" << std::flush;

            MPI_Send(&sleep_times[next_st_idx], 1, MPI_INT, w_rank, WORKTAG, MPI_COMM_WORLD);
            ++next_st_idx;                    // give worker more sleep lol
        } else {
            std::cout << "[root] send STOP to w=" << w_rank << "\n" << std::flush;
            MPI_Send(nullptr, 0, MPI_INT, w_rank, STOPTAG, MPI_COMM_WORLD);
            --sent_count;                         // left over sent count reduces

    }
}
std::cout << "[root] all workers stopped. total_tasks=" << st_size << "\n" << std::flush;
}
