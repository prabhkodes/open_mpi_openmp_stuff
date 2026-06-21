#include "../include/CMatrix.hpp"
#include "../../mpi_timer/include/mpi_timer.hpp"

int main(int argc, char** argv ) {
    MPI_Init(&argc, &argv);

    int  world_size, world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size < 2) {
        std::cerr << "World size must be at least two for this example\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::cout << "Running matrix overload vectors MPI" << std::endl;

    CMatrix<double, 100> A{};
    CMatrix<double, 100> B{};

    A.fill();
    B.fill();

    if (world_rank %2 == 0) {
        {   // timed scope
            CSimple_Timer t_sleep{"Printing from world rank%2=0"};
            write_to_file(A, "A_rank" + std::to_string(world_rank) + ".txt");
        }
    } else {
        {   // timed scope
            CSimple_Timer t_sleep{"Printing from world rank%2=0"};
            write_to_file(B, "B_rank" + std::to_string(world_rank) + ".txt");
        }
    }
    MPI_Finalize();
    return 0;
}
