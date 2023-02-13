#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

enum { SIZE_MATRIX = 16,
    ARBITRARY_VALUE = 0,
    ZERO_VALUE = 0
};

const double τ =  1e-2;
const double ε = 1e-7;

typedef double* dynamicArray;
typedef double* dynamicMatrix;

#define MPI_Allgather MPI_Allgather(multiplyPartMatrix, SIZE_MATRIX/cntProcess, MPI_DOUBLE, vectorUtility, SIZE_MATRIX/cntProcess, MPI_DOUBLE, MPI_COMM_WORLD);

__attribute__((unused)) void PrintMatrix(dynamicArray matrix) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        for(int j = 0; j < SIZE_MATRIX; ++j) {
            std::cout << matrix[i* SIZE_MATRIX + j];
        }
        std::cout << std::endl;
    }
}

void PrintVector(dynamicArray vector, int size) {
    for(int j = 0; j < size; ++j) {
        std::cout << (double)vector[j] << " ";
    }
    printf("\n");
}

dynamicArray GenerateDynamicArray(int size) {
    dynamicMatrix vector = new double[size];
    assert(vector != NULL);
    return vector;
}

dynamicArray GenerateSolutionVector() {
    dynamicArray partMatrix =  GenerateDynamicArray(SIZE_MATRIX);
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        partMatrix[i] = ARBITRARY_VALUE;
    }
    return partMatrix;
}


dynamicArray GenerateVectorRightParts() {
    dynamicArray vector = GenerateDynamicArray(SIZE_MATRIX);
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vector[i] = SIZE_MATRIX+1;
    }
    return vector;
}

void GenerateVectorArbitraryValue(dynamicArray vector, int size, const double value) {
    for(int i = 0; i < size; ++i) {
        vector[i] = value;
    }
}

dynamicMatrix GeneratePartMatrix(const int& rank, const int& countProcess) {
    int partSizeMatrix = SIZE_MATRIX * (SIZE_MATRIX / countProcess);
    dynamicMatrix partMatrix = GenerateDynamicArray(partSizeMatrix);
    GenerateVectorArbitraryValue(partMatrix, partSizeMatrix, 1.0);

    for(int  i = 0, offset = 0; i < SIZE_MATRIX / countProcess; ++i) {
        partMatrix[offset+rank+i] = 2.0;
        offset += SIZE_MATRIX;
    }
    return partMatrix;
}

dynamicArray MultiplyVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, int cntProcess) {
    GenerateVectorArbitraryValue(result, SIZE_MATRIX/cntProcess, ZERO_VALUE);
    for(int i = 0; i < SIZE_MATRIX / cntProcess; ++i) {
        for  (int j = 0; j < SIZE_MATRIX; ++j) {
            result[i] += vector1[j+i*SIZE_MATRIX] * vector2[j];
        }
    }
    return result;
}

dynamicArray MinusVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

dynamicArray MultiplyVectorByConstant(dynamicArray vector, double constant) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vector[i] *= constant;
    }
    return vector;
}

double FormingFirstNorm(const dynamicArray vector) {
    double sumVector = 0;
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        sumVector += vector[i]*vector[i];
    }
    return sqrt(sumVector);
}

double NormCalculation(const dynamicArray vector1, const dynamicArray vector2, dynamicArray vectorUtility) {
    return FormingFirstNorm(MinusVectors(vector1, vector2, vectorUtility)) / FormingFirstNorm(vector2);
}

bool IsFirstNormMoreEpsilon(const dynamicArray vector1, const dynamicArray vector2, dynamicArray vectorUtility) {
    return !(NormCalculation(vector1, vector2, vectorUtility) < ε);
}

void CopyVector(dynamicArray copyVector, const dynamicArray sourceVector) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        copyVector [i] = sourceVector[i];
    }
}

void  DeleteVectors(dynamicMatrix v1, dynamicArray v2, dynamicArray v3, dynamicArray v4, dynamicArray v5, dynamicArray v6) {
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    delete[] v5;
    delete[] v6;
}

//x^(n+1) = x^n – τ(Ax^n – b)
double* IterativeMethod(const int& rank, const int& cntProcess) {
    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess);
    dynamicArray b = GenerateVectorRightParts();
    dynamicArray x = GenerateSolutionVector();

    dynamicMatrix vectorResult = GenerateDynamicArray(SIZE_MATRIX);
    dynamicMatrix multiplyPartMatrix = GenerateDynamicArray(SIZE_MATRIX / cntProcess);
    dynamicMatrix vectorUtility = GenerateDynamicArray(SIZE_MATRIX);
    dynamicMatrix multiplyVectors = GenerateDynamicArray(SIZE_MATRIX);
    do {
        multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, cntProcess);

//        MPI_Allgather(multiplyPartMatrix,
//                      SIZE_MATRIX/cntProcess,
//                      MPI_DOUBLE,
//                      vectorUtility,
//                      SIZE_MATRIX/cntProcess,
//                      MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather;
        CopyVector(multiplyVectors, vectorUtility);
        vectorResult = MinusVectors(x, MultiplyVectorByConstant(MinusVectors(vectorUtility, b, vectorResult), τ), vectorResult); // vectorUtility
        CopyVector(x, vectorResult);
    } while(IsFirstNormMoreEpsilon(multiplyVectors, b, vectorUtility));

    DeleteVectors(A, b, x, multiplyPartMatrix, vectorUtility, multiplyVectors);
    return vectorResult;
}


int main(int argc, char* argv[]) {
    //if(argc != 1) return 1;

    MPI_Init(&argc, &argv);
    int rank = 0, cntProcess = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    double* vector = IterativeMethod(rank, cntProcess);

    MPI_Finalize();

    if(rank == 0)
        PrintVector(vector, SIZE_MATRIX);

    return 0;
}