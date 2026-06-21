#pragma once
#include "../../mpi_methods/include/get_mpi_datatype.hpp"
#include <mpi.h>
#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <thread>
#include <vector>

template <typename T>
void do_send_recv(int world_rank, int world_size, std::size_t n, bool per_element) {
    if (world_size < 2) { return; }

    std::vector<T> recv_buff(n);
    std::vector<T> send_buff;

    MPI_Datatype mpi_dt = get_mpi_datatype<T>();

    if (world_rank == 0) {
        send_buff.resize(n);


        std::cout << "world_rank = " << world_rank
                  << " Sending: " << send_buff.size()
                  << (per_element ? " (element-by-element)\n" : " (bulk)\n");

        if (per_element) {
            for (std::size_t i = 0; i < n; ++i) {
                MPI_Send(send_buff.data() + i, 1, mpi_dt, 1, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Send(send_buff.data(), static_cast<int>(n), mpi_dt, 1, 0, MPI_COMM_WORLD);
        }

    } else { // world_rank == 1
        std::cout << "world_rank = " << world_rank
                  << " Rec: " << recv_buff.size()
                  << (per_element ? " (element-by-element)\n" : " (bulk)\n");

        if (per_element) {
            for (std::size_t i = 0; i < n; ++i) {
                MPI_Recv(recv_buff.data() + i, 1, mpi_dt, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            MPI_Recv(recv_buff.data(), static_cast<int>(n), mpi_dt, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}



