#define Grid double*

enum AREA_CHANGE_SYMBOL {
    Dx = 2,
    Dy = 2,
    Dz = 2,
    x_0 = -1,
    y_0 = -1,
    z_0 = -1
};

enum SIZE_GRID {
    Nx = 532,
    Ny = 230,
    Nz = 204
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
    static constexpr double Hx = ValueStep(AREA_CHANGE_SYMBOL::Dx, SIZE_GRID::Nx);
    static constexpr double Hy =  ValueStep(AREA_CHANGE_SYMBOL::Dy, SIZE_GRID::Ny);
    static constexpr double Hz = ValueStep(AREA_CHANGE_SYMBOL::Dz, SIZE_GRID::Nz);
};

constexpr double CalcCoefficientH(const double H) {
    return 2 / H;
}

struct Const {
    static constexpr double a = 1e5;
    static constexpr double EPLSILOND = 1e-8;
    static constexpr int INITIAL_APPROXIMATION = 0;
    static constexpr int ROOT = 0;
    static constexpr int MPI_TAG_LOW = 24;
    static constexpr int MPI_TAG_UPPER = 42;

    static constexpr double SQUARE_STEP_Hx = GridSteps::Hx * GridSteps::Hx;
    static constexpr double SQUARE_STEP_Hy =  GridSteps::Hy * GridSteps::Hy;
    static constexpr double SQUARE_STEP_Hz = GridSteps::Hz * GridSteps::Hz;

    static constexpr double STEP_Hx = CalcCoefficientH(SQUARE_STEP_Hx);
    static constexpr double STEP_Hy = CalcCoefficientH(SQUARE_STEP_Hy);
    static constexpr double STEP_Hz = CalcCoefficientH(SQUARE_STEP_Hz);

    static constexpr double COEFFICIENT =
            1 / ( STEP_Hx + STEP_Hy + STEP_Hz + Const::a);

    static constexpr int INCREASE_RANK = 1;
};