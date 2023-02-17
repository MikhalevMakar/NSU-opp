#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

enum { SIZE_MATRIX = 1024,
    ARBITRARY_VALUE = 0,
    ZERO_VALUE = 0
};
//τ
//ε
const double τ =  1e-5;
const double ε = 1e-5;

typedef double* dynamicArray;
typedef double* dynamicMatrix;

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

void GenerateVectorArbitraryValue(double* vector, int sizeNew, double value) {
    for(int i = 0; i < sizeNew; ++i) {
        vector[i] = value;
    }
}

dynamicArray GenerateSolutionVector(int fictitiousSize, int sizeOrigin) {
    dynamicArray vector =  GenerateDynamicArray(fictitiousSize);
    GenerateVectorArbitraryValue(vector, fictitiousSize, ZERO_VALUE);
    for(int i = 0; i < sizeOrigin; ++i) {
        vector[i] = ARBITRARY_VALUE;
    }
    return vector;
}

dynamicArray GenerateVectorRightParts(int fictitiousSize, int sizeOrigin) {
    dynamicArray vector = GenerateDynamicArray(fictitiousSize);
    GenerateVectorArbitraryValue(vector, fictitiousSize, ZERO_VALUE);
    for(int i = 0; i < sizeOrigin; ++i) {
        vector[i] = sizeOrigin+1;
    }
    return vector;
}

int GetBalanceSizeVector(const int& sizeOrigin, const int& countProcess) {
    int buildNewSize = sizeOrigin;
    while(buildNewSize % countProcess != 0) {
        ++buildNewSize;
    }
    return buildNewSize;
}

int GetCntCurrentFillLineMatrix(const int& rank, const int& cntProcess, const int& sizeNew) {
    int resultFillLine = 0;
    for (int startRank = 1; startRank < rank+1; ++startRank) {
        resultFillLine += sizeNew / cntProcess;
    }
    return resultFillLine;
}

dynamicMatrix GeneratePartMatrix(const int& rank, const int& cntProcess, const int& fictitiousSize) {
    int partFictitiousSizeMatrix = fictitiousSize * (fictitiousSize / cntProcess);

    dynamicMatrix partMatrix = GenerateDynamicArray(partFictitiousSizeMatrix);

    GenerateVectorArbitraryValue(partMatrix, partFictitiousSizeMatrix, ZERO_VALUE);

    int countRows = fictitiousSize / cntProcess;
    int index = countRows * rank;
    int  numberCntLine = GetCntCurrentFillLineMatrix(rank, cntProcess,  fictitiousSize);

    for (int i = 0, offset = 0; i < countRows; ++i, offset += fictitiousSize) {
        if(numberCntLine+i < SIZE_MATRIX) {
            for (int j = 0; j < SIZE_MATRIX; ++j) {
                partMatrix[offset + j] = 1.0;
            }
            partMatrix[index + i] = 2.0;
        }
        index += fictitiousSize;
    }
    return partMatrix;
}

dynamicArray MultiplyVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, const int& size, const int& cntProcess) {
    GenerateVectorArbitraryValue(result, size/cntProcess, ZERO_VALUE);
    for(int i = 0; i < size / cntProcess; ++i) {
        for  (int j = 0; j < size; ++j) {
            result[i] += vector1[j+i*size] * vector2[j];
        }
    }
    return result;
}

dynamicArray MinusVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, const int& size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

dynamicArray MultiplyVectorByConstant(const dynamicArray vector, const double& constant, dynamicArray result, const int& size) {
    for(int i = 0; i < size; ++i) {
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

bool IsFirstNormMoreEpsilon(const double& v1, const double& v2) {
    return !((v1/v2) < ε);
}

void CopyVector(dynamicArray copyVector, const dynamicArray sourceVector, int size) {
    for(int i = 0; i < size; ++i) {
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
    int fictitiousSize = GetBalanceSizeVector(SIZE_MATRIX, cntProcess);
    int countRows = fictitiousSize / cntProcess;
    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess, fictitiousSize);

    dynamicArray b = GenerateVectorRightParts(fictitiousSize, SIZE_MATRIX);
    dynamicArray x = GenerateSolutionVector(fictitiousSize, SIZE_MATRIX);


    dynamicMatrix vectorResult = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix multiplyPartMatrix = GenerateDynamicArray(countRows);
    dynamicMatrix vectorUtility = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix compareVectors = GenerateDynamicArray(fictitiousSize);

    double normB = FormingFirstNorm(b);

//    do {
//        multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, cntProcess);
//
//        MPI_Allgather(multiplyPartMatrix, fictitiousSize / cntProcess, MPI_DOUBLE, vectorUtility,
//                      fictitiousSize / cntProcess, MPI_DOUBLE, MPI_COMM_WORLD);
//
//        compareVectors = MinusVectors(vectorUtility, b, compareVectors, fictitiousSize);
//
//        vectorResult = MinusVectors(x,
//                                    MultiplyVectorByConstant(compareVectors, τ, vectorResult, fictitiousSize),
//                                    vectorResult, fictitiousSize);
//
//        CopyVector(x, vectorResult, fictitiousSize);
//
//    } while(IsFirstNormMoreEpsilon(FormingFirstNorm(compareVectors), normB));
    do {
        multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, cntProcess);

        MPI_Allgather(multiplyPartMatrix,
                      fictitiousSize / cntProcess,
                      MPI_DOUBLE,
                      vectorUtility,
                      fictitiousSize / cntProcess,
                      MPI_DOUBLE,
                      MPI_COMM_WORLD);

        CopyVector(compareVectors, vectorUtility, fictitiousSize);

        vectorResult = MinusVectors(vectorUtility, b, vectorResult, fictitiousSize);

        vectorResult = MultiplyVectorByConstant(vectorResult, τ, vectorResult, fictitiousSize);

        vectorResult = MinusVectors(x, vectorResult, vectorResult, fictitiousSize);

        CopyVector(x, vectorResult, fictitiousSize);

    } while(IsFirstNormMoreEpsilon(multiplyVectors, b, vectorUtility, fictitiousSize));


    DeleteVectors(A, b, x, multiplyPartMatrix, vectorUtility, compareVectors);
    return vectorResult;
}

int main(int argc, char* argv[]) {
    double startTime, endTime;
    startTime = MPI_Wtime();
    MPI_Init(&argc, &argv);
    int rank = 0, cntProcess = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    dynamicArray vector = IterativeMethod(rank, cntProcess);

    MPI_Finalize();
    endTime = MPI_Wtime();

    if(rank == 0) {
        printf("RESULT VECTOR\n");
        PrintVector(vector, SIZE_MATRIX);
        std::cout << "\nTIME: " << endTime - startTime << std::endl;
    }
    delete[] vector;
    return 0;
}