#include <iostream>
#include <cassert>
#include <mpi.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include "Context.h"

template <typename T, typename... Args>
void OverrideFree(T matrix, Args... matrices) {
    free(matrix);
    if constexpr(sizeof...(matrices) > 0) {
        OverrideFree(matrices...);
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
    double x = GetCoordinate(AREA_CHANGE_Ω::x_0, GridSteps::Hx, i);
    double y = GetCoordinate(AREA_CHANGE_Ω::y_0, GridSteps::Hy, j);
    double z = GetCoordinate(AREA_CHANGE_Ω::z_0, GridSteps::Hz, k);
    return x*x + y*y + z*z;
}

double RightHandSideEquation(const int i, const int j, const int k) {
    return 6 - Const::a * CalculateFunctionValue(i, j, k);
}

Grid GenerateGrid(const std::vector<int> vectorOffset, const std::vector<int> vectorCountLine,
                  const int rank) {
    Grid grid = MemoryAllocatedGrid(vectorCountLine[rank] * SIZE_GRID::Ny * SIZE_GRID::Nz);

    for(int i = 0, beginPosI = vectorOffset[rank]; i < vectorCountLine[rank]; ++i, ++beginPosI)
        for(int j = 0; j < Ny; ++j)
            for(int k = 0; k < Nz; ++k) {
                grid[i * Ny * Nz + j * Nz + k] =
                        (IsZero(beginPosI, j, k) || beginPosI == Nx-1 || j == Ny-1 || k == Nz-1)
                        ?
                        CalculateFunctionValue(beginPosI, j, k) : Const::INITIAL_APPROXIMATION;
            }
    return grid;
}

int GetIndexGrid(const int i, const int j, const int k) {
    return i * SIZE_GRID::Ny * SIZE_GRID::Nz + j * SIZE_GRID::Nz + k;
}

void CalcNextIteration(Grid grid, const Index* index, const int offsetX,
                       const double Fx, const double Fy, const double Fz, double& maximumChange) {
    const double currValue = grid[GetIndexGrid(index->i, index->j, index->k)];

    grid[GetIndexGrid(index->i, index->j, index->k)] =
            (Fx + Fy + Fz - RightHandSideEquation(index->i + offsetX, index->j, index->k)) * Const::COEFFICIENT;

    maximumChange = std::max(maximumChange, fabs(grid[GetIndexGrid(index->i, index->j, index->k)] - currValue));
}

void CalculationInnerPartGrid(Grid grid, double& maximumChange,
                              const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                              const int rank) {
    double Fx, Fy, Fz;

    Index* index = new Index();

    for(int i = 1; i < vectorCountLine[rank]-1; ++i) {
        for(int j = 1; j < Ny-1; ++j) {
            for(int k = 1; k < Nz-1; ++k) {
                *index = {i, j, k};

                Fx = (grid[GetIndexGrid(i + 1, j, k)] + grid[GetIndexGrid(i - 1, j, k)]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[GetIndexGrid(i, j + 1, k)] + grid[GetIndexGrid(i, j - 1, k)]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[GetIndexGrid(i, j, k + 1)] + grid[GetIndexGrid(i, j, k - 1)]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange);
            }
        }
    }
    OverrideFree(index);
}

void CalculationEdgesPartGrid(Grid grid, const Grid bufferLow, const Grid bufferUpper,
                              const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                              double& maximumChange, const int rank, const int countThread) {
    double Fx, Fy, Fz;
    const int indexUpper = 0, indexLower = vectorCountLine[rank] - 1;

    Index* index = new Index();

    for(int j = 1; j < Ny-1; ++j) {
        for (int k = 1; k < Nz - 1; ++k) {
            if (rank != Const::ROOT) {
                *index = {indexUpper, j, k};

                Fx = (grid[GetIndexGrid(indexUpper + 1, j, k)] + bufferLow[j * Nz + k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[GetIndexGrid(indexUpper, j + 1, k)] + grid[GetIndexGrid(indexUpper, j - 1, k)]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[GetIndexGrid(indexUpper, j, k + 1)] + grid[GetIndexGrid(indexUpper, j, k - 1)]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange);
            }

            if (rank != countThread - Const::INCREASE_RANK) {
                *index = {indexLower, j, k};

                Fx = (grid[GetIndexGrid(indexLower - 1, j, k)] + bufferUpper[j * Nz + k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[GetIndexGrid(indexLower, j + 1, k)] + grid[GetIndexGrid(indexLower, j - 1, k)]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[GetIndexGrid(indexLower, j, k + 1)] + grid[ GetIndexGrid(indexLower, j, k - 1)]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange);
            }
        }
    }
    OverrideFree(index);
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

double FindMaxDelta(const Grid grid, const int rank,
                    const std::vector<int>& vectorCountLine, const std::vector<int>& vectorOffset) {
    double maxChange = DBL_MIN, currentCalc;

    for (int i = 0; i < vectorCountLine[rank]; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                currentCalc = CalculateFunctionValue(i + vectorOffset[rank], j, k);
                maxChange = std::max(maxChange, fabs(grid[GetIndexGrid(i, j, k)] - currentCalc));
            }
        }
    }
    return maxChange;
}

void ICommutation(const int rank, const int countThread, const int offsetX, Grid bufferReceivedLow,
                  Grid bufferReceivedUpper, Grid grid, MPI_Request* request)  {
    if(rank != Const::ROOT) {

        MPI_Isend(grid, SIZE_GRID::Ny * SIZE_GRID::Nz, MPI_DOUBLE,
                  rank - Const::INCREASE_RANK, Const::MPI_TAG_UPPER, MPI_COMM_WORLD, &request[0]);

        MPI_Irecv(bufferReceivedLow, SIZE_GRID::Ny * SIZE_GRID::Nz, MPI_DOUBLE,
                  rank - Const::INCREASE_RANK, Const::MPI_TAG_LOW, MPI_COMM_WORLD, &request[1]);
    }

    if( rank != countThread-1) {

        MPI_Isend(grid + offsetX * SIZE_GRID::Ny * SIZE_GRID::Nz, SIZE_GRID::Ny * SIZE_GRID::Nz, MPI_DOUBLE,
                  rank + Const::INCREASE_RANK, Const::MPI_TAG_LOW, MPI_COMM_WORLD, &request[2]);

        MPI_Irecv(bufferReceivedUpper, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank + Const::INCREASE_RANK, Const::MPI_TAG_UPPER, MPI_COMM_WORLD, &request[3]);
    }
}

void WaitingRequest(MPI_Request* request, const int rank, const int countThread) {
    if (rank != 0) {
        MPI_Wait(&request[0], MPI_STATUS_IGNORE);
        MPI_Wait(&request[1], MPI_STATUS_IGNORE);
    }

    if(rank != countThread-1) {
        MPI_Wait(&request[2], MPI_STATUS_IGNORE);
        MPI_Wait(&request[3], MPI_STATUS_IGNORE);
   }
}

void FindMaxChange(double* vectorMaxChange, double& maxChange, const int countThread) {
    MPI_Allgather(&maxChange, 1, MPI_DOUBLE, vectorMaxChange, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    for(int i = 0; i < countThread; ++i) {
        maxChange = std::max(maxChange, vectorMaxChange[i]);
    }
}

void CalculateMethodJacobi(Grid grid, const int rank, const int countThread,
                           const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset) {
    MPI_Request requestSendRecv[4];
    double maximumChange = DBL_MIN;

    double* vectorMaxChange = MemoryAllocatedGrid(countThread);

    Grid bufferReceivedLow = MemoryAllocatedGrid(SIZE_GRID::Ny * SIZE_GRID::Nz);
    Grid bufferReceivedUpper = MemoryAllocatedGrid(SIZE_GRID::Ny * SIZE_GRID::Nz);

    do {
        maximumChange = DBL_MIN;
        ICommutation(rank, countThread, vectorCountLine[countThread-1]-1,
                     bufferReceivedLow, bufferReceivedUpper, grid, requestSendRecv);

        CalculationInnerPartGrid(grid, maximumChange, vectorCountLine, vectorOffset, rank);

        WaitingRequest(requestSendRecv, rank, countThread);

        CalculationEdgesPartGrid(grid, bufferReceivedLow, bufferReceivedUpper,
                                 vectorCountLine, vectorOffset,
                                 maximumChange, rank, countThread);

        FindMaxChange(vectorMaxChange, maximumChange, countThread);
    } while(maximumChange >= Const::ε);

    maximumChange = FindMaxDelta(grid, rank, vectorCountLine, vectorOffset);
    FindMaxChange(vectorMaxChange, maximumChange, countThread);

    if(rank == Const::ROOT)
        std::cout << "CHANGE_VALUE: " << maximumChange << std::endl;

    OverrideFree(bufferReceivedUpper, bufferReceivedLow, vectorMaxChange);
}

void RunMethodJacobi() {
    int rank, countThread;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    std::vector<int> vectorCountLine, vectorOffset;

    GenerateVectorOffset(vectorCountLine, vectorOffset, countThread);
    Grid grid = GenerateGrid(vectorOffset, vectorCountLine, rank);

    double startTime = MPI_Wtime();
    CalculateMethodJacobi(grid, rank, countThread, vectorCountLine, vectorOffset);
    double endTime = MPI_Wtime();

    if(rank == Const::ROOT)
        std::cout << std::endl << "TIME: " << endTime - startTime << " seconds" << std::endl;

    OverrideFree(grid);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    RunMethodJacobi();

    MPI_Finalize();
    return 0;
}