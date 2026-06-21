#include <iostream>
#include <cstdlib> 
#include <vector>
#include <mpi.h>


template <typename T, int N>
class CMatrix{
    std::vector<T> data;
    int world_size{};
    int world_rank{};
    MPI_Comm comm{MPI_COMM_NULL};
public:
    CMatrix();
    
    void fill();

    void print() const;

    template<typename M, int K>
    friend std::ostream& operator<< (std::ostream& outputstream, const CMatrix<M, K>& matrix);

};

// Constructor of CMatrix<T, N>
template<typename T, int N>
CMatrix<T, N>::CMatrix() : data(N * N, T{}), 
            world_size(MPI_Comm_size(MPI_COMM_WORLD, &world_size)), 
            world_rank(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank)), 
            comm(MPI_COMM_WORLD)  {} 


// Fill method of CMatrix<T, N>
template <typename T, int N>
void CMatrix<T, N>::fill() {
        for (int i=0; i<N; i++) 
            for (int j=0; j<N; j++)
               data[i * N + j] = static_cast<T>(rand()) / RAND_MAX * 10.0;
    }


// Print method of CMatrix<T, N>
template <typename T, int N>
void CMatrix<T,N>::print() const {
        std::cout << "Im printing from: " << world_rank << std::endl;
        for(int i=0; i<N; i++){
            std::cout << "| ";
            for (int j=0; j<N; j++)
                std::cout << data[i * N + j] << " ";
            std::cout << "| " << std::endl;
        }
    }


// OPERATOR OVERLOADERS
template <typename M, int K>
std::ostream& operator<< (std::ostream& outputstream, const CMatrix<M, K>& matrix) {
        std::cout << "Im printing from: " << matrix.world_rank << std::endl;
        for (int i=0; i<K; i++) {
            outputstream << "| ";
            for (int j=0; j<K; j++)
                outputstream << matrix.data[i * K + j] << " ";
            outputstream << "|" << std::endl;
        }
        return outputstream;
    }

#include <fstream>
#include <stdexcept>

template <typename M, int K>
void write_to_file(const CMatrix<M, K>& matrix, const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Failed to open file: " + filename);
    ofs << matrix;  
    ofs.close();
}

