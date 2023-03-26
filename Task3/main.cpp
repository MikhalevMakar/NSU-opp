#include <iostream>#include <mpi.h>#define Matrix double*enum COORDS {    X,    Y};enum PROCESS {    p1 = 4,    p2 = 2};enum SIZE_MATRIX {    n1 = 15,    n2 = 10,    n3 = 20,};enum {    NUMBER_DIMS = 2,    ROOT = 0,    ZERO_VALUE = 0};void PrintMatrix(Matrix matrix, int sizeColumn, int sizeRow) {      for (int i = 0; i < sizeRow; ++i) {          for (int j = 0; j < sizeColumn; ++j) {              std::cout << matrix[i*sizeColumn+j] << " ";          }          std::cout << std::endl;      }}Matrix GenerateMatrix(int sizeColumn, int sizeRow) {    Matrix matrix = new double[sizeColumn * sizeRow];    for (int i = 0; i < sizeRow; ++i) {        for (int j = 0; j < sizeColumn; ++j) {            matrix[i * sizeColumn + j] = static_cast<double>(i * sizeColumn + j);        }    }    return matrix;}void CreationCommunicators(int size, int rank, MPI_Comm& comm2d, MPI_Comm& commColumn, MPI_Comm& commRow) {    int dims[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    MPI_Dims_create(size, NUMBER_DIMS, dims);    int periods[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    int reorder = 1;    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_DIMS, dims, periods, reorder, &comm2d);    int coords[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    MPI_Cart_coords(comm2d, rank, NUMBER_DIMS, coords);    MPI_Comm_split(comm2d, coords[COORDS::X], coords[COORDS::Y], &commColumn);    MPI_Comm_split(comm2d, coords[COORDS::Y], coords[COORDS::X], &commRow);}bool CheckResult(Matrix result, Matrix rightDecision, int column, int row) {    for(int i = 0; i < row; ++i) {        for(int j = 0; j < column; ++j) {            if(result[i*column+j] != rightDecision[i*column+j])                return false;        }    }    return true;}void RunMultiplicationMatrix(int rank, int cntProcess) {    Matrix A; Matrix B; Matrix C;    if(rank == ROOT) {        A = GenerateMatrix(SIZE_MATRIX::n1, SIZE_MATRIX::n2);        B = GenerateMatrix(SIZE_MATRIX::n2, SIZE_MATRIX::n3);        C = new double[SIZE_MATRIX::n1*SIZE_MATRIX::n3];    }    Matrix partMatrixA = new double[SIZE_MATRIX::n1*SIZE_MATRIX::n2/p1];    Matrix partMatrixB = new double[SIZE_MATRIX::n2/p2*SIZE_MATRIX::n3];    Matrix partMatrixC = new double[SIZE_MATRIX::n2/p2 * SIZE_MATRIX::n3/p1];    MPI_Comm comm2d, commColumn, commRow;    CreationCommunicators(cntProcess, rank, comm2d, commColumn, commRow);}void FreeMatrix(Matrix m1, Matrix m2, Matrix m3) {    free(m1);    free(m2);    free(m3);}void FreeComm(MPI_Comm* comm1, MPI_Comm* comm2, MPI_Comm* comm3) {    MPI_Comm_free(comm1);    MPI_Comm_free(comm2);    MPI_Comm_free(comm3);}int main(int argc, char* argv[]) {    MPI_Init(&argc, &argv);    double startTime = MPI_Wtime();    int rank = 0, cntProcess = 0;    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);    RunMultiplicationMatrix(rank, cntProcess);    double endTime = MPI_Wtime();    if(rank == ROOT) {        std::cout << "\nTIME: " << endTime - startTime << std::endl;    }    MPI_Finalize();    return EXIT_SUCCESS;}