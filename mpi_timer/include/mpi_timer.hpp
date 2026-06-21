#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <cstdint>
#include <vector>
#include <mpi.h>
#include <fstream>

struct TimeTable {
    std::string   label;
    std::uint64_t time_elapsed;
};

class CSimple_Timer {
    std::string what;
    std::chrono::steady_clock::time_point start{};
    int world_rank = 0, world_size = 1;

    inline static std::vector<TimeTable> table; // shared per process

public:
    explicit CSimple_Timer(const std::string& label) : what(label) {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        start = std::chrono::steady_clock::now();
    }
    ~CSimple_Timer() {
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::steady_clock::now() - start);
        table.push_back({what, static_cast<std::uint64_t>(us.count())});
    }

static void gather_timings(const std::string& path = "timings.txt") {
    int root{0}, rank{0}, size{1};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    // convert the TimeTable to a string buffer for writing later
    std::string file_o_string;
    file_o_string.reserve(table.size() * 48);
    for (const auto& row : table) {
        file_o_string += row.label;
        file_o_string += '\t';
        file_o_string += std::to_string(row.time_elapsed);
        file_o_string += '\n';
    }

    // Bytes to send would be string size
    const int send_n = static_cast<int>(file_o_string.size());

    // Gather recv_buffer and buffers on root
    std::vector<int> recv_buffer(rank == root ? size : 0);
    MPI_Gather(&send_n, 1, MPI_INT,
               (rank == root ? recv_buffer.data() : nullptr), 1, MPI_INT,
               root, MPI_COMM_WORLD);
    

    // Compute displacements and total on root
    std::vector<int> displs;
    int total = 0;
    if (rank == root) {
        displs.resize(size);
        int off = 0;
        for (int i = 0; i < size; ++i) { displs[i] = off; off += recv_buffer[i]; }
        total = off;
    }

    // Gather payloads (bytes) on root
    std::vector<char> all(rank == root ? total : 0);
    MPI_Gatherv(file_o_string.data(), send_n, MPI_CHAR,
                (rank == root ? all.data()   : nullptr),
                (rank == root ? recv_buffer.data() : nullptr),
                (rank == root ? displs.data(): nullptr),
                MPI_CHAR, root, MPI_COMM_WORLD);


    // Root writes grouped by rank, preserving lines
    if (rank == root) {
        std::ofstream out(path, std::ios::app);
        if (out) {
            for (int r = 0; r < size; ++r) {
                out << "rank " << r << ":\n";
                if (recv_buffer[r] > 0) {
                    const char* p = all.data() + displs[r];
                    out.write(p, recv_buffer[r]);
                    if (p[recv_buffer[r] - 1] != '\n') out << '\n';
                }
            }
            out.flush();
        }
    }
}

};
