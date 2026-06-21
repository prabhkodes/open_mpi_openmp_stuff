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
void do_scatter(int world_rank, int world_size,
                    std::size_t n_elems, bool as_bytes) {

    std::vector<T> recv_buff(n_elems);
    std::vector<T> send_buff;

    const MPI_Datatype dtype = as_bytes ? MPI_BYTE : get_mpi_datatype<T>();
    
    const std::size_t count_sz = as_bytes ? n_elems * sizeof(T) : n_elems;

    const int count = static_cast<int>(count_sz); 

    if (world_rank == 0) {

        send_buff.resize(n_elems * static_cast<std::size_t>(world_size));
        std::cout << "world_rank = " << world_rank
                  << (as_bytes ? " Sending (bytes): " : " Sending (elems): ")
                  << (as_bytes ? send_buff.size() * sizeof(T) : send_buff.size())
                  << std::endl;
    } else { // world_rank == 1
        std::cout << "world_rank = " << world_rank
                  << (as_bytes ? " Rec (bytes): " : " Rec (elems): ")
                  << (as_bytes ? recv_buff.size() * sizeof(T) : recv_buff.size())
                  << std::endl;
    }

    MPI_Scatter(
        world_rank == 0 ? static_cast<void*>(send_buff.data()) : nullptr,
        count, dtype,
        static_cast<void*>(recv_buff.data()),
        count, dtype,
        0, MPI_COMM_WORLD
    );
}
