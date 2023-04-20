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
    Ny = 4,
    Nz = 4
};

struct Index {
    int i;
    int j;
    int k;
};

constexpr double ValueStep(const int D, const int N) {
    return static_cast<double>(D / static_cast<double>(N - 1));
}

struct GridSteps {
    static constexpr double Hx = ValueStep(AREA_Ω::Dx, SIZE_GRID::Nx);
    static constexpr double Hy =  ValueStep(AREA_Ω::Dy, SIZE_GRID::Ny);
    static constexpr double Hz = ValueStep(AREA_Ω::Dz, SIZE_GRID::Nz);
};

constexpr double CalcCoefficientH(const double H) {
    return 2 / (H*H);
}

struct Const {
    static constexpr double a = 1e5;
    static constexpr double ε = 1e-8;
    static constexpr int INITIAL_APPROXIMATION = 0;
    static constexpr int ROOT = 0;
    static constexpr int MPI_TAG_LOW = 24;
    static constexpr int MPI_TAG_UPPER = 42;

    static constexpr double STEP_Hx = CalcCoefficientH(GridSteps::Hx);
    static constexpr double STEP_Hy = CalcCoefficientH(GridSteps::Hy);
    static constexpr double STEP_Hz = CalcCoefficientH(GridSteps::Hz);

    static constexpr double SQUARE_STEP_Hx = GridSteps::Hx * GridSteps::Hx;
    static constexpr double SQUARE_STEP_Hy =  GridSteps::Hy * GridSteps::Hy;
    static constexpr double SQUARE_STEP_Hz = GridSteps::Hz * GridSteps::Hz;

    static constexpr double COEFFICIENT =
            1 / (STEP_Hx + STEP_Hy + STEP_Hz + Const::a);

};

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
    double x = GetCoordinate(AREA_Ω::x_0, GridSteps::Hx, i);
    double y = GetCoordinate(AREA_Ω::y_0, GridSteps::Hy, j);
    double z = GetCoordinate(AREA_Ω::z_0, GridSteps::Hz, k);
    return x*x + y*y + z*z;
}

double RightHandSideEquation(const int i, const int j, const int k) {
    return 6 - Const::a * CalculateFunctionValue(i, j, k);
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
                        CalculateFunctionValue(beginPosI, j, k) : Const::INITIAL_APPROXIMATION;
            }
    return grid;
}

void CalcNextIteration(Grid grid, const Index* index, const int offsetX,
                       const int Fx, const int Fy, const int Fz, double& maximumChange) {

    const double currValue = grid[index->i * Ny * Nz + index->j * Nz + index->k];

    grid[index->i * Ny * Nz + index->j * Nz + index->k] =
            (Fx + Fy + Fz - RightHandSideEquation(index->i + offsetX, index->j, index->k)) * Const::COEFFICIENT;

    if (maximumChange > fabs(grid[index->i * Ny * Nz + index->j * Nz + index->k] - currValue))
        maximumChange = fabs(grid[index->i * Ny * Nz + index->j * Nz + index->k] - currValue);
}

void CalculationInnerPartGrid(Grid grid, double& maximumChange,
                              const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                              const int rank) {
//    double currValue, Fx, Fy, Fz;
//
//    for(int i = 1; i < vectorCountLine[rank]-1; ++i) {
//        for (int j = 1; j < Ny - 1; ++j) {
//            for (int k = 1; k < Nz - 1; ++k) {
//
//                Fx = (grid[(i + 1) * Ny * Nz + j * Nz + k] + grid[(i - 1) * Ny * Nz + j * Nz + k]) / Const::SQUARE_STEP_Hx;
//                Fy = (grid[i * Ny * Nz + (j + 1) * Nz + k] + grid[i * Ny * Nz + (j - 1) * Nz + k]) / Const::SQUARE_STEP_Hy;
//                Fz = (grid[i * Ny * Nz + j * Nz + k + 1] + grid[i * Ny * Nz + j * Nz + k - 1]) / Const::SQUARE_STEP_Hz;
//
//                currValue = grid[i * Ny * Nz + j * Nz + k];
//
//                grid[i * Ny * Nz + j * Nz + k] =
//                        (Fx + Fy + Fz - RightHandSideEquation(i+vectorOffset[rank], j, k)) * Const::COEFFICIENT;
//
//                if (maximumChange > fabs(grid[i * Ny * Nz + j * Nz + k] - currValue))
//                    maximumChange = fabs(grid[i * Ny * Nz + j * Nz + k] - currValue);
//            }
//        }
//    }

    double currValue, Fx, Fy, Fz;

    for(int i = 1; i < Nx-1; ++i) {
        for(int j = 1; j < Ny-1; ++j) {
            for(int k = 1; k < Nz-1; ++k) {

                Fx = (grid[(i+1) * Ny * Nz + j * Nz + k] + grid[(i-1) * Ny * Nz + j * Nz + k]) / (GridSteps::Hx*GridSteps::Hx);
                Fy = (grid[i * Ny * Nz + (j+1) * Nz + k] + grid[i * Ny * Nz + (j-1) * Nz + k]) / (GridSteps::Hy*GridSteps::Hy);
                Fz = (grid[i * Ny * Nz + j * Nz + k+1] + grid[i * Ny * Nz + j * Nz + k-1]) / (GridSteps::Hz*GridSteps::Hz);

                currValue = grid[i * Ny * Nz + j * Nz + k];

                grid[i * Ny * Nz + j * Nz + k] = (Fx + Fy + Fz - RightHandSideEquation(i+vectorOffset[rank], j, k)) * Const::COEFFICIENT;
                if(maximumChange > fabs(grid[i * Ny * Nz + j * Nz + k] - currValue))
                    maximumChange = fabs(grid[i * Ny * Nz + j * Nz + k] - currValue);
            }
        }
    }
}

void CalculationEdgesPartGrid(Grid grid, const Grid bufferLow, const Grid bufferUpper,
                              const std::vector<int> vectorCountLine, const std::vector<int> vectorOffset,
                              double& maximumChange, const int rank, const int countThread) {
    double Fx, Fy, Fz;
    const int indexUpper = 0, indexLower = vectorCountLine[rank]-1;

    Index* index = new Index();

    for(int j = 1; j < Ny-1; ++j)
        for(int k = 1; k < Nz-1; ++k) {
            if (rank != Const::ROOT) {
                *index = {indexUpper, j, k};

                Fx = (grid[indexUpper*Ny*Nz + Nz*j + k] + bufferLow[j*Nz+k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[indexUpper*Ny*Nz + (j+1) * Nz + k] + grid[indexUpper*Ny*Nz + (j-1)*Nz + k]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[indexUpper*Ny*Nz + j * Nz + k + 1] + grid[indexUpper*Ny*Nz + j*Nz + k - 1]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange);
            }

            if(rank != countThread-1) {
                *index = {indexLower, j, k};

                Fx = (grid[indexLower*Ny*Nz + Nz*j + k] + bufferUpper[j*Nz + k]) / Const::SQUARE_STEP_Hx;
                Fy = (grid[indexLower*Ny*Nz + (j+1) * Nz + k] + grid[indexLower*Ny*Nz + (j-1) * Nz + k]) / Const::SQUARE_STEP_Hy;
                Fz = (grid[indexLower*Ny*Nz + j*Nz + k + 1] + grid[indexLower*Ny*Nz + j*Nz + k - 1]) / Const::SQUARE_STEP_Hz;

                CalcNextIteration(grid, index, vectorOffset[rank],
                                  Fx, Fy, Fz, maximumChange);
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

double FindMaxDelta(const Grid grid, std::vector<int>& vectorOffset, const int rank) {
    double maxChange = DBL_MIN, currentCalc;

    for (int i = 1; i < Nx-1; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            for (int k = 1; k < Nz-1; ++k) {
                currentCalc = CalculateFunctionValue(i+vectorOffset[rank], j, k);

                if(maxChange < fabs(grid[i * Ny * Nz + j * Nz + k] - currentCalc)) {
                    maxChange = fabs(grid[i * Ny * Nz + j * Nz + k] - currentCalc);
                }
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
            std::cout << buffer[i*Nz+j] << " ";
        }
        printf("\n");
    }
    printf("\n\n");
}

void RunMethodJacobi(const int rank) {
    int countThread;
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    double maximumChange = DBL_MAX;
    std::vector<int> vectorCountLine, vectorOffset;

    Grid bufferReceivedLow = MemoryAllocatedGrid(SIZE_GRID::Ny*SIZE_GRID::Nz);
    Grid bufferReceivedUpper = MemoryAllocatedGrid(SIZE_GRID::Ny*SIZE_GRID::Nz);

    GenerateVectorOffset(vectorCountLine, vectorOffset, countThread);
    Grid grid = GenerateGrid(vectorOffset, vectorCountLine, rank, countThread);

    MPI_Request request[4];

    do {
        ICommutation(rank, countThread, vectorCountLine[countThread-1]-1,
                     bufferReceivedLow, bufferReceivedUpper, grid, request);

        CalculationInnerPartGrid(grid, maximumChange, vectorCountLine, vectorOffset, rank);

        WaitingRequest(request, rank, countThread);

        CalculationEdgesPartGrid(grid, bufferReceivedLow, bufferReceivedUpper,
                                 vectorCountLine, vectorOffset,
                                 maximumChange, rank, countThread);

    } while(maximumChange > Const::ε);

    MPI_Barrier(MPI_COMM_WORLD);

    PrintGrid(grid, vectorCountLine, rank);

    if(rank == Const::ROOT)
        std::cout << "CHANGE_VALUE: " << FindMaxDelta(grid, vectorOffset, rank) << std::endl;

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