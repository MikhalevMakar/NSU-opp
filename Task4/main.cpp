#include <iostream>
#include <assert.h>

enum AREA_Ω{
    Dx = 4,
    Dy = 4,
    Dz = 4,
    x_0 = 4,
    y_0 = 4,
    z_0 = 4
};

enum SIZE_GRID {
    Nx = 5,
    Ny = 5,
    Nz = 5
};

constexpr int ValueStep(int D, int N) {
    return D / (N - 1);
}

enum GRID_STEPS {
    Hx = ValueStep(AREA_Ω::Dx, SIZE_GRID::Nx),
    Hy =  ValueStep(AREA_Ω::Dy, SIZE_GRID::Ny),
    Hz = ValueStep(AREA_Ω::Dz, SIZE_GRID::Nz)
};

#define Grid double*

Grid MemoryAllocatedGrid() {
    Grid grid = new double[SIZE_GRID::Nx *
                           SIZE_GRID::Ny *
                           SIZE_GRID::Nz];
    assert(grid);
    return grid;
}


template <typename T, typename... Args>
bool IsZero(T value, Args... digits) {
    if(value == 0) return true;
    if constexpr(sizeof...(digits) > 0) {
        return IsZero(digits...);
    }
    return false;
}



int FirstNorma(int x, int y, int z) {
    return x*x + y*y + z*z;
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

Grid GenerateGrid() {
    Grid grid = MemoryAllocatedGrid();

   for(int x = 0; x < Nx; ++x) {
    for(int y = 0; y < Ny; ++y) {
        for(int z = 0; z < Nz; ++z) {
            if(x == 0 || y == 0 || z == 0 || x == Nx-1 || y == Ny-1 || z == Nz-1) {
                grid[x * Ny * Nz + y * Nz + z] = FirstNorma(x, y, z);
            } else {
                grid[x * Ny * Nz + y * Nz + z] = 0;
            }
        }
    }

}
    return grid;
}


int main() {

    Grid grid = GenerateGrid();
    PrintGrid(grid);

    return 0;
}