
#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <cfloat>
#include <cmath>

#define Grid double*
#define Vector double*

enum AREA_Ω{
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

constexpr double ValueStep(int D, int N) {
    return D / (double)(N - 1);
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

    constexpr double COEFFICIENT_Hx = CalcCoefficientH(GRID_STEPS::Hx);
    constexpr double COEFFICIENT_Hy = CalcCoefficientH(GRID_STEPS::Hy);
    constexpr double COEFFICIENT_Hz = CalcCoefficientH(GRID_STEPS::Hz);
}

void PrintGrid(Grid grid) {
    for(int x = 0; x < Nx; ++x) {
        for(int y = 0; y < Ny; ++y) {
            for(int z = 0; z < Nz; ++z) {
                std::cout << grid[x * Ny * Nz + y * Nz + z] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
}

Grid MemoryAllocatedGrid() {
    Grid grid = new double[SIZE_GRID::Nx *
                           SIZE_GRID::Ny *
                           SIZE_GRID::Nz];
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

Grid GenerateGrid() {
    Grid grid = MemoryAllocatedGrid();

    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                grid[i * Ny * Nz + j * Nz + k] =
                        (IsZero(i, j, k) || i == Nx-1 || j == Ny-1 || k == Nz-1)
                        ?
                        CalculateFunctionValue(i, j, k) : CONST::INITIAL_APPROXIMATION;
            }
        }
    }
    return grid;
}

void CalculationRequiredFunction(Grid grid, double& maximumChange) {
    double currValue, Fx, Fy, Fz;

    for(int i = 1; i < Nx-1; ++i) {
        for(int j = 1; j < Ny-1; ++j) {
            for(int k = 1; k < Nz-1; ++k) {

                Fx = (grid[(i+1) * Ny * Nz + j * Nz + k] + grid[(i-1) * Ny * Nz + j * Nz + k]) / (GRID_STEPS::Hx*GRID_STEPS::Hx);
                Fy = (grid[i * Ny * Nz + (j+1) * Nz + k] + grid[i * Ny * Nz + (j-1) * Nz + k]) / (GRID_STEPS::Hy*GRID_STEPS::Hy);
                Fz = (grid[i * Ny * Nz + j * Nz + k+1] + grid[i * Ny * Nz + j * Nz + k-1]) / (GRID_STEPS::Hz*GRID_STEPS::Hz);

                currValue = grid[i * Ny * Nz + j * Nz + k];

                double coefficient = (1 / (CONST::COEFFICIENT_Hx + CONST::COEFFICIENT_Hy + CONST::COEFFICIENT_Hz + CONST::a));

                grid[i * Ny * Nz + j * Nz + k] = (Fx + Fy + Fz - RightHandSideEquation(i, j, k)) * coefficient;
                maximumChange = std::max(maximumChange, fabs(grid[i * Ny * Nz + j * Nz + k] - currValue));
                if(maximumChange > fabs(grid[i * Ny * Nz + j * Nz + k] - currValue))

                    maximumChange = fabs(grid[i * Ny * Nz + j * Nz + k] - currValue);
            }
        }
    }
};

void LaunchFindGridDecisions(Grid grid) {
    // PrintGrid(grid);
    double maximumChange = DBL_MAX;
    int cnt = 0;
    do {
    CalculationRequiredFunction(grid, maximumChange);
    //} while(maximumChange > CONST::ε);
    } while(++cnt < 1);
}

double FindMaxChange(Grid grid) {
    double maxChange = DBL_MIN;
    for (int i = 0; i < Nx; ++i) {
       for (int j = 0; j < Ny; ++j) {
           for (int k = 0; k < Nz; ++k) {
                if(maxChange < fabs(grid[i * Ny * Nz + j * Nz + k]-CalculateFunctionValue(i, j, k)))
                    maxChange = fabs(grid[i * Ny * Nz + j * Nz + k]-CalculateFunctionValue(i, j, k));
            }
        }
    }
    return maxChange;
}



void RunMethodJacobi() {
    int cntProcess, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Grid grid = NULL;

    if(rank == CONST::ROOT)  grid = GenerateGrid();

    std::cout << "CHANGE_VALUE: " << FindMaxChange(grid) << std::endl;

    double startTime = MPI_Wtime();

    LaunchFindGridDecisions(grid);

    double endTime = MPI_Wtime();


    PrintGrid(grid);
    std::cout << "CHANGE_VALUE: " << FindMaxChange(grid) << std::endl;
    std::cout << std::endl << "TIME: " << endTime - startTime << " seconds"<< std::endl;
    free(grid);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    RunMethodJacobi();

    MPI_Finalize();
    return 0;
}