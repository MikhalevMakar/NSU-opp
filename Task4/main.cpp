#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <cfloat>
#include <cmath>
#include <vector>

#define Grid double*
enum AREA_Ω {
    Dx = 2,
    Dy = 2,
    Dz = 2,
    x_0 = -1,
    y_0 = -1,
    z_0 = -1
};

enum SIZE_GRID {
    Nx = 4,
    Ny = 5,
    Nz = 5
};

constexpr double ValueStep(const int D, const int N) {
    return static_cast<double>(D / static_cast<double>(N - 1));
}

namespace GRID_STEPS {
    constexpr double Hx = ValueStep(AREA_Ω::Dx, SIZE_GRID::Nx);
    constexpr double Hy =  ValueStep(AREA_Ω::Dy, SIZE_GRID::Ny);
    constexpr double Hz = ValueStep(AREA_Ω::Dz, SIZE_GRID::Nz);
};

constexpr double CalcCoefficientH(double H) {
    return 2 / (H*H);
}

namespace CONST {
    constexpr double a = 1e5;
    constexpr double ε = 1e-8;
    constexpr int INITIAL_APPROXIMATION = 0;
    constexpr int ROOT = 0;
    constexpr int MPI_TAG = 1234;

    constexpr double COEFFICIENT_Hx = CalcCoefficientH(GRID_STEPS::Hx);
    constexpr double COEFFICIENT_Hy = CalcCoefficientH(GRID_STEPS::Hy);
    constexpr double COEFFICIENT_Hz = CalcCoefficientH(GRID_STEPS::Hz);

    constexpr double COEFFICIENT =
            1 / (CONST::COEFFICIENT_Hx + CONST::COEFFICIENT_Hy + CONST::COEFFICIENT_Hz + CONST::a);
}

void PrintGrid(Grid grid, const std::vector<int>& vectorCountLine, const int rank) {
    for(int x = 0; x < vectorCountLine[rank]; ++x) {
        for(int y = 0; y < Ny; ++y) {
            for(int z = 0; z < Nz; ++z) {
                std::cout << grid[x * Ny * Nz + y * Nz + z] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
}

Grid MemoryAllocatedGrid(int size) {
    Grid grid = new double[size];
    assert(grid);
    return grid;
}

template <typename T, typename... Args>
bool IsZero(const T value, const Args... digits) {
    if(value == 0) return true;
    if constexpr(sizeof...(digits) > 0) {
        return IsZero(digits...);
    }
    return false;
}

double GetCoordinate(const double startCoord, const double step, const int index) {
    return startCoord + index*step;
}

double CalculateFunctionValue(const int i, const int j, const int k) {
    double x = GetCoordinate(AREA_Ω::x_0, GRID_STEPS::Hx, i);
    double y = GetCoordinate(AREA_Ω::y_0, GRID_STEPS::Hy, j);
    double z = GetCoordinate(AREA_Ω::z_0, GRID_STEPS::Hz, k);
    return x*x + y*y + z*z;
}

double RightHandSideEquation(const int i, const int j, const int k) {
    return 6 - CONST::a * CalculateFunctionValue(i, j, k);
}

Grid GenerateGrid(const std::vector<int> vectorOffset, const std::vector<int> vectorCountLine,
                  const int rank, const int countThread) {
    Grid grid = MemoryAllocatedGrid(SIZE_GRID::Nx * SIZE_GRID::Ny * SIZE_GRID::Nz / countThread);

    for(int i = 0, beginPosI = vectorOffset[rank]; i < vectorCountLine[rank]; ++i, ++beginPosI)
        for(int j = 0; j < Ny; ++j)
            for(int k = 0; k < Nz; ++k) {
                grid[i * Ny * Nz + j * Nz + k] =
                        (IsZero(beginPosI, j, k) || beginPosI == Nx-1 || j == Ny-1 || k == Nz-1)
                        ?
                        CalculateFunctionValue(beginPosI, j, k) : CONST::INITIAL_APPROXIMATION;
            }
    return grid;
}

void CalculationInnerPartGrid(Grid grid, double& maximumChange,
                                 const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                                 const int rank) {
    double currValue, Fx, Fy, Fz;

    for(int i = 1; i < vectorCountLine[rank]-1; ++i) {
        for(int j = 1; j < Ny-1; ++j) {
            for(int k = 1; k < Nz-1; ++k) {

                Fx = (grid[(i+1) * Ny * Nz + j * Nz + k] + grid[(i-1) * Ny * Nz + j * Nz + k]) / (GRID_STEPS::Hx*GRID_STEPS::Hx);
                Fy = (grid[i * Ny * Nz + (j+1) * Nz + k] + grid[i * Ny * Nz + (j-1) * Nz + k]) / (GRID_STEPS::Hy*GRID_STEPS::Hy);
                Fz = (grid[i * Ny * Nz + j * Nz + k+1] + grid[i * Ny * Nz + j * Nz + k-1]) / (GRID_STEPS::Hz*GRID_STEPS::Hz);

                currValue = grid[i * Ny * Nz + j * Nz + k];

                grid[i * Ny * Nz + j * Nz + k] = (Fx + Fy + Fz - RightHandSideEquation(i+vectorOffset[rank], j, k)) * CONST::COEFFICIENT;
                if(maximumChange > fabs(grid[i * Ny * Nz + j * Nz + k] - currValue))
                    maximumChange = fabs(grid[i * Ny * Nz + j * Nz + k] - currValue);
            }
        }
    }
}

void GenerateVectorOffset(std::vector<int>& vectorCountLine, std::vector<int>& vectorOffset, const int size) {
    vectorCountLine.resize(size, SIZE_GRID::Nx / size);
    vectorOffset.resize(size);

    int remainder = SIZE_GRID::Nx % size;
    for(int i = 0; i < remainder; ++i) {
        ++vectorCountLine[i];
    }
    for (int i = 0, offset = 0; i < size; offset += vectorCountLine[i++]) {
        vectorOffset[i] = offset;
    }
}

double FindMaxDelta(Grid grid) {
    double maxChange = DBL_MIN;
    for (int i = 0; i < Nx; ++i)
       for (int j = 0; j < Ny; ++j)
           for (int k = 0; k < Nz; ++k) {
                if(maxChange < fabs(grid[i * Ny * Nz + j * Nz + k]-CalculateFunctionValue(i, j, k)))
                    maxChange = fabs(grid[i * Ny * Nz + j * Nz + k]-CalculateFunctionValue(i, j, k));
           }
}

void print(Grid arr) {
    for (int i = 0; i < SIZE_GRID::Ny; ++i) {
        for (int j = 0; j < SIZE_GRID::Nz; ++j) {
            std::cout << arr[i*SIZE_GRID::Ny + j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void ICommutation(const int rank, const int countThread, const int offsetX, Grid bufferReceivedLow, Grid bufferReceivedUpper, Grid grid)  {
    MPI_Request request[4];

    if(rank != CONST::ROOT) {

        MPI_Isend(grid, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank-1, 0, MPI_COMM_WORLD, &request[0]);

        MPI_Isend(bufferReceivedLow, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank-1, 1, MPI_COMM_WORLD, &request[1]);

    } else if( rank != countThread-1) {

        MPI_Isend(grid+offsetX*SIZE_GRID::Ny*SIZE_GRID::Nz, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank+1, 1, MPI_COMM_WORLD, &request[2]);

        MPI_Isend(bufferReceivedUpper, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank+1, 0, MPI_COMM_WORLD, &request[3]);
    }
}

void RunMethodJacobi(const int rank) {
    int countThread;
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    double maximumChange = DBL_MAX;
    std::vector<int> vectorCountLine, vectorOffset;

    Grid bufferReceivedLow = MemoryAllocatedGrid(SIZE_GRID::Ny * SIZE_GRID::Nz);
    Grid bufferReceivedUpper = MemoryAllocatedGrid(SIZE_GRID::Ny * SIZE_GRID::Nz);

    GenerateVectorOffset(vectorCountLine, vectorOffset, countThread);
    Grid grid = GenerateGrid(vectorOffset, vectorCountLine, rank, countThread);

    ICommutation(rank, countThread, vectorOffset[countThread-1], bufferReceivedLow, bufferReceivedUpper, grid);

    printf("\nLOWER\n");
    print(bufferReceivedLow);

    printf("\nUPPER\n");
    print(bufferReceivedUpper);

//    do {
//        CalculationInnerPartGrid(grid, maximumChange, vectorCountLine, vectorOffset, rank);
//    } while(maximumChange > CONST::ε);

    MPI_Barrier(MPI_COMM_WORLD);
    PrintGrid(grid, vectorCountLine, rank);

    if(rank == CONST::ROOT)
        std::cout << "CHANGE_VALUE: " << FindMaxDelta(grid) << std::endl;

     free(grid);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double startTime = MPI_Wtime();

    RunMethodJacobi(rank);

    double endTime = MPI_Wtime();

    if(rank == CONST::ROOT)
        std::cout << std::endl << "TIME: " << endTime - startTime << " seconds"<< std::endl;

    MPI_Finalize();
    return 0;
}

//#include <mpi.h>
//#include <iostream>
//
//int main(int argc, char** argv) {
//    MPI_Init(&argc, &argv);
//
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    int buffer[10];
//
//    if (rank == 0) {
//        // Заполняем буфер данными
//        for (int i = 0; i < 10; i++) {
//            buffer[i] = i;
//        }
//
//        // Отправляем данные процессу с рангом 1
//        MPI_Request request;
//        MPI_Isend(buffer, 10, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
//
//        // Ждем завершения передачи данных
//        MPI_Wait(&request, MPI_STATUS_IGNORE);
//    } else if (rank == 1) {
//        // Получаем данные от процесса с рангом 0
//        MPI_Request request;
//        MPI_Irecv(buffer, 10, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
//
//        // Ждем завершения приема данных
//        MPI_Wait(&request, MPI_STATUS_IGNORE);
//
//        // Выводим полученные данные
//        for (int i = 0; i < 10; i++) {
//            std::cout << buffer[i] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//    MPI_Finalize();
//    return 0;
//}
