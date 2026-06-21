#include "../include/include.hpp"
#include <mpi.h>
#include <thread>
#include <chrono>

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Test the timer by thread sleep op
    {CSimple_Timer t_sleep{"Sleep 50000 micro s"}; std::this_thread::sleep_for(std::chrono::microseconds(50000)); }

    // 2 MiB for testing
    std::size_t n_ints    = (2ull * 1024 * 1024) / sizeof(int);
    std::size_t n_doubles = (2ull * 1024 * 1024) / sizeof(double);

    { CSimple_Timer t{"Scatter<int> full vector (native)"}; do_scatter<int>(world_rank, world_size, n_ints, false); }
    { CSimple_Timer t{"SendRecv<int> full vector"};   do_send_recv<int>(world_rank, world_size, n_ints, false); }
    { CSimple_Timer t{"SendRecv<int> by element"};  do_send_recv<int>(world_rank, world_size, n_ints, true); }
    { CSimple_Timer t{"Scatter<int> as bytes"};   do_scatter<int>(world_rank, world_size, n_ints, true); }

    { CSimple_Timer t{"Scatter<double> full vector (native)"}; do_scatter<double>(world_rank, world_size, n_doubles, false); }
    { CSimple_Timer t{"SendRecv<double> full vector"}; do_send_recv<double>(world_rank, world_size, n_doubles, false); }
    { CSimple_Timer t{"SendRecv<double> by element"};  do_send_recv<double>(world_rank, world_size, n_doubles, true); }
    { CSimple_Timer t{"Scatter<double> as bytes"};  do_scatter<double>(world_rank, world_size, n_doubles, true); }

    { CSimple_Timer t{"PUT<double> by element"};  do_put<double>(world_rank, world_size, n_doubles, true); }
    { CSimple_Timer t{"PUT<double> full vector"}; do_put<double>(world_rank, world_size, n_doubles, false); }

    CSimple_Timer::gather_timings("timings.txt");

    MPI_Finalize();
    return 0;
}
