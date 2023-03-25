#include <iostream>
#include <mpi.h>

void CreationCommunicators() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int numberDims = 2;

    int dims[2] = {0, 0};
    MPI_Dims_create(size, numberDims, dims);
    int sizeX = dims[0];
    int sizeY = dims[1];
    std::cout << sizeX << " " << sizeY << std::endl;

    MPI_Comm comm2d;
    int periods[2] = {0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, numberDims, dims, periods, reorder, &comm2d);

    int rank;
    MPI_Comm_rank(comm2d, &rank);

    int coords[2];
    MPI_Cart_get(comm2d, 2, dims, periods, coords);
    int rankY = coords[0];
    int rankX = coords[1];
    std::cout << rankX << " " << rankY << std::endl;

}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    double startTime = MPI_Wtime();

    int rank = 0, cntProcess = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    CreationCommunicators();

    double endTime = MPI_Wtime();

    if(rank == 0) {
        std::cout << "\nTIME: " << endTime - startTime << std::endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}