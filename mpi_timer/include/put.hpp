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
void do_put(int world_rank, int world_size, std::size_t n, bool per_element) {
    if (world_size < 2) { return; }

    std::vector<T> recv_buff(n);   
    std::vector<T> send_buff;      

    MPI_Datatype mpi_dt = get_mpi_datatype<T>();

    
    MPI_Win win;
    MPI_Win_create(
        recv_buff.data(),                   
        n * sizeof(T),                      // size in bytes
        static_cast<int>(sizeof(T)),        // disp_unit
        MPI_INFO_NULL,
        MPI_COMM_WORLD,
        &win
    );

    if (world_rank == 0) {
        send_buff.resize(n);
        

        std::cout << "world_rank = " << world_rank
                  << (per_element ? " PUT element-by-element\n" : " PUT bulk\n");

        
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 1, 0, win);

        if (per_element) {
            for (std::size_t i = 0; i < n; ++i) {
                
                MPI_Put(send_buff.data() + i, 1, mpi_dt,
                        1, static_cast<MPI_Aint>(i), 1, mpi_dt, win);
            }
        } else {
            
            MPI_Put(send_buff.data(), static_cast<int>(n), mpi_dt,
                    1, 0, static_cast<int>(n), mpi_dt, win);
        }

        
        MPI_Win_unlock(1, win);

    } else { // world_rank == 1 (target)
        std::cout << "world_rank = " << world_rank
                  << " waiting for PUT into recv_buff of size " << recv_buff.size() << "\n";
    }



    MPI_Win_free(&win);
}