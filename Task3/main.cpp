#include <iostream>
#include <mpi.h>

#define Matrix double*

enum COORDS {
    X,
    Y
};

enum SIZE_MATRIX {
    n1 = 8,
    n2 = 6,
    n3 = 12
};

enum {
    NUMBER_DIMS = 2,
    ROOT = 0
};

template <typename T, typename... Args>
void FreeMatrices(T matrix, Args... matrices) {
    free(matrix);
    if constexpr(sizeof...(matrices) > 0) {
        FreeMatrices(matrices...);
    }
}

template <typename T, typename... Args>
void FreeCommunicators(T comm, Args... communicators) {
    MPI_Comm_free(&comm);
    if constexpr(sizeof...(communicators) > 0) {
        FreeCommunicators(communicators...);
    }
}

template <typename T, typename... Args>
void  FreeDatatypeMPI(T type, Args... types) {
      MPI_Type_free(&type);
     if constexpr(sizeof...(types) > 0) {
        FreeDatatypeMPI(types...);
    }
}

void PrintMatrix(const Matrix matrix, const int sizeRow, const int sizeColumn) {
      for (int i = 0; i < sizeRow; ++i) {
          for (int j = 0; j < sizeColumn; ++j) {
              std::cout << matrix[i*sizeColumn+j] << " ";
          }
          std::cout << std::endl;
      }
}

Matrix GenerateMatrix(const int sizeRow, const int sizeColumn) {
    Matrix matrix = new double[sizeColumn * sizeRow];
    for (int i = 0; i < sizeRow; ++i) {
        for (int j = 0; j < sizeColumn; ++j) {
            matrix[i * sizeColumn + j] = static_cast<double>(i * sizeColumn + j + 1);
        }
    }
    return matrix;
}

void CreationCommunicators(const int size, const int rank, int coords[], int dims[],
                           MPI_Comm* comm2d, MPI_Comm* commColumn, MPI_Comm* commRow) {
    MPI_Dims_create(size, NUMBER_DIMS, dims);

    int periods[NUMBER_DIMS] = {0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_DIMS, dims, periods, reorder, comm2d);
    MPI_Cart_coords(*comm2d, rank, NUMBER_DIMS, coords);
    int coordsRow[] = {0, 1};
    int coordsColumn[] = {1, 0};
    MPI_Cart_sub(*comm2d, coordsRow, commRow);
    MPI_Cart_sub(*comm2d, coordsColumn, commColumn);
}

void CreateParametersGatherv(int recvCounts[], int displs[], int dims[]) {
    int sizeShift = SIZE_MATRIX::n1  * dims[COORDS::Y] / dims[COORDS::X];
     for (int i = 0; i < dims[COORDS::X]; ++i) {
         for (int j = 0; j < dims[COORDS::Y]; ++j) {
             recvCounts[i * dims[COORDS::Y] + j] = 1;
             displs[i * dims[COORDS::Y] + j] =  i * sizeShift  + j;
         }
     }
}


void CreateMatrixCFromSeparatedBlock(const Matrix partC, Matrix C,
                                     const int sizePartRow, const int sizePartColumn,
                                     const int countColumn, int dims[],
                                     const int rank, const int cntProcess,
                                     const MPI_Comm comm2d) {
    MPI_Datatype vectorType, newType;
    MPI_Type_vector(sizePartRow, sizePartColumn, countColumn, MPI_DOUBLE, &vectorType);
    MPI_Type_commit(&vectorType);

    MPI_Type_create_resized(vectorType, 0, sizePartColumn * sizeof(double), &newType);
    MPI_Type_commit(&newType);

    int receiveCounts[cntProcess];
    int displays[cntProcess];

    if(rank == ROOT)
        CreateParametersGatherv(receiveCounts, displays, dims);

    MPI_Gatherv(partC, sizePartColumn * sizePartRow, MPI_DOUBLE, C, receiveCounts, displays, newType, ROOT, comm2d);

    FreeDatatypeMPI(vectorType, newType);
}

void MatrixMultiplication(const Matrix matrixFirst, const Matrix matrixSecond, Matrix result,
                          const int N, const int M, const int K) {
    for(int i = 0; i < N*K; ++i) {
        result[i] = 0;
    }
    for(int i = 0; i < N; ++i) {
        for(int k = 0; k < M; ++k) {
            for(int j = 0; j < K; ++j) {
                result[i * K + j] += matrixFirst[i * M + k] * matrixSecond[k * K + j];
            }
        }
    }
}


bool CheckResult(const Matrix result, const Matrix rightDecision,
                 const int row, const int column) {
    for(int i = 0; i < row; ++i) {
        for(int j = 0; j < column; ++j) {
            if(result[i*column+j] != rightDecision[i*column+j])
                return false;
        }
    }
    return true;
}

void DistributionMatrixIntoNodes(const Matrix matrixA, const Matrix matrixB,
                                 Matrix partMatrixA, Matrix partMatrixB,
                                 const MPI_Comm commColumn, const MPI_Comm commRow,
                                 int dims[], int coords[]) {

    if(coords[COORDS::Y] == 0) {
        MPI_Scatter(matrixA, (SIZE_MATRIX::n1 * SIZE_MATRIX::n2) / dims[COORDS::X], MPI_DOUBLE, partMatrixA,
                    (SIZE_MATRIX::n1 * SIZE_MATRIX::n2) / dims[COORDS::X], MPI_DOUBLE, ROOT, commColumn);
    }

    if(coords[COORDS::X] == 0) {
        MPI_Datatype vectorType, newType;

        MPI_Type_vector(SIZE_MATRIX::n2, SIZE_MATRIX::n3 / dims[COORDS::Y], SIZE_MATRIX::n3, MPI_DOUBLE, &vectorType);
        MPI_Type_commit(&vectorType);

        MPI_Type_create_resized(vectorType, 0, SIZE_MATRIX::n3 / dims[COORDS::Y] * sizeof(double), &newType);
        MPI_Type_commit(&newType);

        MPI_Scatter(matrixB, 1, newType, partMatrixB,
                    (SIZE_MATRIX::n3 * SIZE_MATRIX::n2) / dims[COORDS::Y], MPI_DOUBLE, ROOT, commRow);
        FreeDatatypeMPI(vectorType, newType);
    }

    MPI_Bcast(partMatrixA, (SIZE_MATRIX::n1 * SIZE_MATRIX::n2) / dims[COORDS::X], MPI_DOUBLE, ROOT, commRow);
    MPI_Bcast(partMatrixB, (SIZE_MATRIX::n3 * SIZE_MATRIX::n2) / dims[COORDS::Y], MPI_DOUBLE, ROOT, commColumn);
}

void LaunchFastMultiplicationMatrices(Matrix A, Matrix B, Matrix C, const int rank) {
    int cntProcess = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    MPI_Comm comm2d, commColumn, commRow;
    int coords[NUMBER_DIMS] = {0, 0};
    int dims[NUMBER_DIMS] = {0, 0};
    CreationCommunicators(cntProcess, rank, coords, dims, &comm2d, &commColumn, &commRow);


    Matrix partMatrixA = new double[(SIZE_MATRIX::n1 * SIZE_MATRIX::n2) / dims[COORDS::X]];
    Matrix partMatrixB = new double[(SIZE_MATRIX::n2 * SIZE_MATRIX::n3) / dims[COORDS::Y]];
    Matrix partMatrixC = new double[(SIZE_MATRIX::n2 * SIZE_MATRIX::n3) / (dims[COORDS::X] * dims[COORDS::Y])];

    DistributionMatrixIntoNodes(A, B,
                                partMatrixA, partMatrixB,
                                commColumn, commRow,
                                dims, coords);

    MatrixMultiplication(partMatrixA, partMatrixB, partMatrixC,
                         SIZE_MATRIX::n1 / dims[COORDS::X],
                         SIZE_MATRIX::n2,
                         SIZE_MATRIX::n3 / dims[COORDS::Y]);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 1)
        PrintMatrix(partMatrixC, n1/dims[X], n3/dims[Y]);


    CreateMatrixCFromSeparatedBlock(partMatrixC,
                                    C,
                                    SIZE_MATRIX::n1 / dims[COORDS::X],
                                    SIZE_MATRIX::n3 / dims[COORDS::Y],
                                    SIZE_MATRIX::n3, dims,
                                    rank, cntProcess,
                                    comm2d);

    FreeMatrices(partMatrixA, partMatrixB, partMatrixC);
    FreeCommunicators(comm2d, commColumn, commRow);
}

void RunMultiplicationMatrix() {
    Matrix A = NULL; Matrix B = NULL; Matrix C = NULL; Matrix CorrectResult = NULL;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT) {
        A = GenerateMatrix(SIZE_MATRIX::n1, SIZE_MATRIX::n2);
        B = GenerateMatrix(SIZE_MATRIX::n2, SIZE_MATRIX::n3);
        C = new double[SIZE_MATRIX::n1 * SIZE_MATRIX::n3];
        CorrectResult = new double[SIZE_MATRIX::n1 * SIZE_MATRIX::n3];
        MatrixMultiplication(A, B, CorrectResult, SIZE_MATRIX::n1, SIZE_MATRIX::n2, SIZE_MATRIX::n3);
    }

    double startTime = MPI_Wtime();

    LaunchFastMultiplicationMatrices(A, B, C, rank);

    double endTime = MPI_Wtime();

    if(rank == ROOT && !CheckResult(C, CorrectResult, SIZE_MATRIX::n1, SIZE_MATRIX::n3))
        std::cout << "Matrix multiplication is wrong\n";
    else if(rank == ROOT)
      PrintMatrix(C, SIZE_MATRIX::n1, SIZE_MATRIX::n3);

    if(rank == ROOT) {
        std::cout << std::endl << "TIME: " << endTime - startTime << std::endl;
        FreeMatrices(A, B, C, CorrectResult);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    RunMultiplicationMatrix();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
