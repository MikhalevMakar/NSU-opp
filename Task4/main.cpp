#include <iostream>
#include <cassert>
#include <mpi.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include "Context.h"

void PrintGrid(const Grid grid, const std::vector<int>& vectorCountLine, const int rank) {
    std::cout << "\n RANK "<< rank << "\n";
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

template <typename T, typename... Args>
void OverrideFree(T matrix, Args... matrices) {
    free(matrix);
    if constexpr(sizeof...(matrices) > 0) {
        OverrideFree(matrices...);
    }
}

double CalculateFunctionValue(const int i, const int j, const int k) {
    double x = GetCoordinate(AREA_CHANGE_SYMBOL::x_0, GridSteps::Hx, i);
    double y = GetCoordinate(AREA_CHANGE_SYMBOL::y_0, GridSteps::Hy, j);
    double z = GetCoordinate(AREA_CHANGE_SYMBOL::z_0, GridSteps::Hz, k);
    return x*x + y*y + z*z;
}

double RightHandSideEquation(const int i, const int j, const int k) {
    return 6 - Const::a * CalculateFunctionValue(i, j, k);
}

Grid GenerateGrid(const std::vector<int> vectorOffset, const std::vector<int> vectorCountLine,
                  const int rank, const int countThread) {
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
                       const double Fx, const double Fy, const double Fz, double& maximumChange, int rank) {
    const double currValue = grid[GetIndexGrid(index->i, index->j, index->k)];

    grid[GetIndexGrid(index->i, index->j, index->k)] =
            (Fx + Fy + Fz - RightHandSideEquation(index->i + offsetX, index->j, index->k)) * Const::COEFFICIENT;

    maximumChange = std::max(maximumChange, fabs(grid[GetIndexGrid(index->i, index->j, index->k)] - currValue));

    if(rank == 0) {
       //std::cout << fabs(grid[index->i * Ny * Nz + index->j * Nz + index->k] - currValue)  << " " << maximumChange << "\n";
    }
}


// .
// .
// .

void CalculationInnerPartGrid(Grid grid, double& maximumChange,
                              const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                              const int rank) {
    double Fx, Fy, Fz;

    Index* index = new Index();

    for(int i = 1; i < vectorCountLine[rank]-1; ++i) {
        for(int j = 1; j < Ny-1; ++j) {
            for(int k = 1; k < Nz-1; ++k) {

                *index = {i, j, k};

                if((i + 1) * Ny * Nz + j * Nz + k > (vectorCountLine[rank]-1) * Ny * Nz + (Ny-1) * Nz + Nz-1) {
                    std::cout << "POSSITION " << (i + 1) * Ny * Nz + j * Nz + k;
                }

                Fx = (grid[(i+1) * Ny * Nz + j * Nz + k] + grid[(i-1) * Ny * Nz + j * Nz + k]) / (GridSteps::Hx*GridSteps::Hx);
                Fy = (grid[i * Ny * Nz + (j+1) * Nz + k] + grid[i * Ny * Nz + (j-1) * Nz + k]) / (GridSteps::Hy*GridSteps::Hy);
                Fz = (grid[i * Ny * Nz + j * Nz + k+1] + grid[i * Ny * Nz + j * Nz + k-1]) / (GridSteps::Hz*GridSteps::Hz);

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange, rank); //delete rank
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

                Fx = (grid[(indexUpper + 1) * Ny * Nz + Nz * j + k] + bufferLow[j * Ny + k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[indexUpper * Ny * Nz + (j + 1) * Nz + k] + grid[indexUpper * Ny * Nz + (j - 1) * Nz + k]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[indexUpper * Ny * Nz + j * Nz + k + 1] + grid[indexUpper * Ny * Nz + j * Nz + k - 1]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange, rank); // delete rank
            }

            if (rank != countThread-1) {
                *index = {indexLower, j, k};

                Fx = (grid[GetIndexGrid(indexLower - 1, j, k)] + bufferUpper[j * Ny + k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[GetIndexGrid(indexLower, j + 1, k)] + grid[GetIndexGrid(indexLower, j - 1, k)]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[GetIndexGrid(indexLower, j, k + 1)] + grid[ GetIndexGrid(indexLower, j, k - 1)]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange, rank); //delete rank
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

double FindMaxDelta(const Grid grid, std::vector<int>& vectorCountLine, std::vector<int>& vectorOffset, const int rank) {
    double maxChange = DBL_MIN, currentCalc;

    for (int i = 0; i < vectorCountLine[rank]; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                currentCalc = CalculateFunctionValue(i+vectorOffset[rank], j, k);
                maxChange = std::max(maxChange, fabs(grid[GetIndexGrid(i, j, k)] - currentCalc));
            }
        }
    }
    return maxChange;
}

void ICommutation(const int rank, const int countThread, const int offsetX, Grid bufferReceivedLow,
                  Grid bufferReceivedUpper, Grid grid, MPI_Request request[])  {
    if(rank != Const::ROOT) {

        MPI_Isend(grid, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank-1, Const::MPI_TAG_UPPER, MPI_COMM_WORLD, &request[0]);

        MPI_Irecv(bufferReceivedLow, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank-1, Const::MPI_TAG_LOW, MPI_COMM_WORLD, &request[1]);
    }

    if( rank != countThread-1) {

        MPI_Isend(grid + offsetX*SIZE_GRID::Ny*SIZE_GRID::Nz, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank+1, Const::MPI_TAG_LOW, MPI_COMM_WORLD, &request[2]);

        MPI_Irecv(bufferReceivedUpper, SIZE_GRID::Ny*SIZE_GRID::Nz, MPI_DOUBLE,
                  rank+1, Const::MPI_TAG_UPPER, MPI_COMM_WORLD, &request[3]);
    }
}

void WaitingRequest(MPI_Request request[], const int rank, const int countThread) {
    if (rank != 0) {
        MPI_Wait(&request[0], MPI_STATUS_IGNORE);
        MPI_Wait(&request[1], MPI_STATUS_IGNORE);
    }

    if(rank != countThread-1) {
        MPI_Wait(&request[2], MPI_STATUS_IGNORE);
        MPI_Wait(&request[3], MPI_STATUS_IGNORE);
   }
}

void print(double* buffer) {
    for(int i = 0; i < Ny; ++i) {
        for(int j = 0; j < Nz; ++j) {
            std::cout << buffer[i*Ny+j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void FindMaxChange(double* vectorMaxChange, double& maxChange, const int countThread) {
    MPI_Allgather(&maxChange, 1, MPI_DOUBLE, vectorMaxChange, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    for(int i = 0; i < countThread; ++i) {
        maxChange = std::max(maxChange, vectorMaxChange[i]);
    }
}

void RunMethodJacobi(const int rank) {
    int countThread;
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    double maximumChange;
    std::vector<int> vectorCountLine, vectorOffset;

    Grid bufferReceivedLow = MemoryAllocatedGrid(SIZE_GRID::Ny*SIZE_GRID::Nz);
    Grid bufferReceivedUpper = MemoryAllocatedGrid(SIZE_GRID::Ny*SIZE_GRID::Nz);

    GenerateVectorOffset(vectorCountLine, vectorOffset, countThread);
    Grid grid = GenerateGrid(vectorOffset, vectorCountLine, rank, countThread);

    MPI_Request request[4];

    double* vectorMaxChange = MemoryAllocatedGrid(countThread);

    do {
        maximumChange = DBL_MIN;
        ICommutation(rank, countThread, vectorCountLine[countThread-1]-1,
                     bufferReceivedLow, bufferReceivedUpper, grid, request);

        CalculationInnerPartGrid(grid, maximumChange, vectorCountLine, vectorOffset, rank);

        WaitingRequest(request, rank, countThread);

        CalculationEdgesPartGrid(grid, bufferReceivedLow, bufferReceivedUpper,
                                 vectorCountLine, vectorOffset,
                                 maximumChange, rank, countThread);

        FindMaxChange(vectorMaxChange, maximumChange, countThread);
    } while(maximumChange >= Const::EPLSILOND);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << rank << " CHANGE_VALUE: " << FindMaxDelta(grid, vectorCountLine, vectorOffset, rank) << std::endl;

    OverrideFree(grid, bufferReceivedUpper, bufferReceivedLow);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double startTime = MPI_Wtime();

    RunMethodJacobi(rank);

    double endTime = MPI_Wtime();

    if(rank == Const::ROOT)
        std::cout << std::endl << "TIME: " << endTime - startTime << " seconds"<< std::endl;

    MPI_Finalize();
    return 0;
}