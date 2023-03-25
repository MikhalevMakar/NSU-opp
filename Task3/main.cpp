#include <iostream>#include <mpi.h>#define Matrix double*enum Coord {    X,    Y};enum {    n1 = 15,    n2 = 15,    n3 = 15,    NUMBER_DIMS = 2,    ROOT = 0,    ZERO_VALUE = 0};Matrix GenerateMatrix(int sizeColumn, int sizeRow) {    Matrix matrix = new double[sizeColumn * sizeRow];    for (int i = 0; i < sizeRow; ++i) {        for (int j = 0; j < sizeColumn; ++j) {            matrix[i * sizeColumn + j] = i * sizeColumn + j;        }    }    return matrix;}void CreationCommunicators(int size, int rank,  MPI_Comm& commColumn, MPI_Comm& commRow) {    int dims[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    MPI_Dims_create(size, NUMBER_DIMS, dims);    MPI_Comm comm2d;    int periods[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    int reorder = 1;    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_DIMS, dims, periods, reorder, &comm2d);    int coords[NUMBER_DIMS] = {ZERO_VALUE, ZERO_VALUE};    MPI_Cart_coords(comm2d, rank, NUMBER_DIMS, coords);    for(int v : coords) {        std::cout << v <<   std::endl;    }    MPI_Comm_split(comm2d, coords[X], coords[Y], &commColumn);    MPI_Comm_split(comm2d, coords[Y], coords[X], &commRow);}bool CheckResult(Matrix result, Matrix rightDecision, int row, int column) {    for(int i = 0; i < row; ++i) {        for(int j = 0; j < column; ++j) {            if(result[i*column+j] != rightDecision[i*column+j])                return false;        }    }    return true;}void FreeMatrix(Matrix m1, Matrix m2, Matrix m3) {    free(m1);    free(m2);    free(m3);}void FreeComm(MPI_Comm* comm1, MPI_Comm* comm2, MPI_Comm* comm3) {    MPI_Comm_free(comm1);    MPI_Comm_free(comm2);    MPI_Comm_free(comm3);}int main(int argc, char* argv[]) {    MPI_Init(&argc, &argv);    double startTime = MPI_Wtime();    int rank = 0, cntProcess = 0;    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);    Matrix C;    if(rank == ROOT) {        Matrix A = GenerateMatrix(n1, n2);        Matrix B = GenerateMatrix(n2, n3);        C = new double[n1*n3];    }    MPI_Comm commColumn, commRow;    CreationCommunicators(cntProcess, rank, commColumn, commRow);    double endTime = MPI_Wtime();    if(rank == ROOT) {        std::cout << "\nTIME: " << endTime - startTime << std::endl;    }    MPI_Finalize();    return EXIT_SUCCESS;}