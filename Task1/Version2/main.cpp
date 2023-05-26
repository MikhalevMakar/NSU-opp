#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

enum { SIZE_VECTOR = 2048,
       ARBITRARY_VALUE = 0,
       SIZE_ONE = 1,
       ZERO_VALUE = 0,
       ROOT = 0
};

const double tau =  1e-5;
const double epsilon = 1e-7;

typedef double* dynamicVector;
typedef double* dynamicMatrix;


void PrintVector(dynamicVector vector, const int size) {
    for(int j = 0; j < size; ++j) {
        std::cout << (double)vector[j] << " ";
    }
    printf("\n");
}

dynamicVector GenerateDynamicArray(const int size) {
    dynamicMatrix vector = new double[size];
    assert(vector != NULL);
    return vector;
}

void GenerateVectorArbitraryValue(dynamicVector vector, const int size, const double value) {
    for(int i = 0; i < size; ++i) {
        vector[i] = value;
    }
}

int GetBalanceSizeVector(const int sizeOrigin, const int countProcess) {
    int buildNewSize = sizeOrigin;
    while(buildNewSize % countProcess != 0) {
        ++buildNewSize;
    }
    return buildNewSize;
}

int GetNumberFillLine(const int rank, const int newSize) {
    int resultFillLine = 0;
    for (int startRank = 1; startRank <= rank; ++startRank) {
        resultFillLine += newSize;
    }
    return resultFillLine;
}


dynamicVector GenerateSolutionVector(const int newSize, const int rank) {
    dynamicVector vector =  GenerateDynamicArray(newSize);

    int  numberCntLine = GetNumberFillLine(rank, newSize);

    for(int i = 0; i < newSize; ++i) {
        vector[i] = (numberCntLine+i < SIZE_VECTOR) ? ARBITRARY_VALUE : ZERO_VALUE;
    }
    return vector;
}

dynamicVector GenerateVectorRightParts(const int newSize, const int rank) {
    dynamicVector vector = GenerateDynamicArray(newSize);

    int  numberCntLine = GetNumberFillLine(rank, newSize);
    for(int i = 0; i < newSize; ++i) {
        vector[i] = (numberCntLine+i < SIZE_VECTOR) ? SIZE_VECTOR + 1 : ZERO_VALUE;
    }
    return vector;
}

int GetCntCurrentFillLineMatrix(const int rank, const int cntProcess, const int newSize) {
    int resultFillLine = 0;
    for (int startRank = 1; startRank <= rank; ++startRank) {
        resultFillLine += newSize / cntProcess;
    }
    return resultFillLine;
}

dynamicMatrix GeneratePartMatrix(const int rank, const int cntProcess, const int fictitiousSize) {
    int partFictitiousSizeMatrix = fictitiousSize * (fictitiousSize / cntProcess);

    dynamicMatrix partMatrix = GenerateDynamicArray(partFictitiousSizeMatrix);

    GenerateVectorArbitraryValue(partMatrix, partFictitiousSizeMatrix, ZERO_VALUE);

    int countRows = fictitiousSize / cntProcess;
    int index = countRows * rank;
    int numberCntLine = GetCntCurrentFillLineMatrix(rank, cntProcess, fictitiousSize);

    for (int i = 0, offset = 0; i < countRows; ++i, offset += fictitiousSize) {
        if (numberCntLine + i + 1 <= SIZE_VECTOR) {
            for (int j = 0; j < SIZE_VECTOR; ++j) {
                partMatrix[offset + j] = 1.0;
            }
            partMatrix[index + i] = 2.0;
        }
        index += fictitiousSize;
    }
    return partMatrix;
}

void MultiplyVectors(const dynamicVector vector1,
                              const dynamicVector vector2,
                              dynamicVector result,
                              const int& size,
                              const int& rank,
                              const int& cntProcess) {
    GenerateVectorArbitraryValue(result, size, ZERO_VALUE);
    for(int i = 0; i < size/cntProcess; ++i) {
        for  (int j = 0; j < size; ++j) {
            result[j] += vector1[j+i*size] * vector2[i];
        }
    }
}

void MinusVectors(const dynamicVector vector1, const dynamicVector vector2,
                           dynamicVector result, const int size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
}

void MultiplyVectorByConstant(const dynamicVector vector, dynamicVector result,
                                       const double constant, const int size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector[i] * constant;
    }
}

double FormingEuclideanNorm(const dynamicVector vector, const int& size) {
    double sumVector = 0;
    for(int i = 0; i < size; ++i) {
        sumVector += vector[i]*vector[i];
    }
    return sumVector;
}

bool IsFirstNormMoreEpsilon(const double v1, const double v2) {
    return !((v1/v2) < epsilon);
}

void CopyVector(dynamicVector copyVector, const dynamicVector sourceVector, const int& size) {
    for(int i = 0; i < size; ++i) {
        copyVector [i] = sourceVector[i];
    }
}

void  DeleteVectors(dynamicVector v1, dynamicVector v2, dynamicVector v3,
                    dynamicVector v4, dynamicVector v5, dynamicVector v6) {
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    delete[] v5;
    delete[] v6;
}

//x^(n+1) = x^n – tau(Ax^n – b)
double* IterativeMethod(const int rank, const int cntProcess) {
    int fictitiousSize = GetBalanceSizeVector(SIZE_VECTOR, cntProcess);
    int fixedSizePartVector = fictitiousSize / cntProcess;

    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess, fictitiousSize);
    dynamicMatrix multiplyMatrixVector = GenerateDynamicArray(fictitiousSize);

    dynamicVector x = GenerateSolutionVector(fixedSizePartVector, rank);

    dynamicVector b = GenerateVectorRightParts(fixedSizePartVector, rank);


    dynamicMatrix vectorResult = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix vectorUtility = GenerateDynamicArray(fixedSizePartVector);

    double normPartB, normB, findPartNorm, resultNorm;
    normPartB = FormingEuclideanNorm(b, fixedSizePartVector);
    MPI_Allreduce(&normPartB, &normB, SIZE_ONE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    dynamicMatrix multiplyPartMatrixVector = GenerateDynamicArray(fictitiousSize);
    do {
        MultiplyVectors(A, x, multiplyPartMatrixVector, fictitiousSize, rank, cntProcess);

        MPI_Allreduce(multiplyPartMatrixVector, multiplyMatrixVector,
                      fictitiousSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        MPI_Scatter(multiplyMatrixVector, fixedSizePartVector, MPI_DOUBLE, vectorUtility,
                    fixedSizePartVector, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

        MinusVectors(vectorUtility, b, vectorUtility, fixedSizePartVector);

        findPartNorm = FormingEuclideanNorm(vectorUtility, fixedSizePartVector);

        MPI_Allreduce(&findPartNorm, &resultNorm, SIZE_ONE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        MultiplyVectorByConstant(vectorUtility, vectorUtility, tau, fixedSizePartVector);

        MinusVectors(x, vectorUtility, vectorUtility, fixedSizePartVector);

        CopyVector(x, vectorUtility, fixedSizePartVector);

    } while(IsFirstNormMoreEpsilon(sqrt(resultNorm), sqrt(normB)));

    MPI_Gather(vectorUtility, fixedSizePartVector, MPI_DOUBLE, vectorResult,
               fixedSizePartVector, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    DeleteVectors(A, b, vectorUtility, multiplyPartMatrixVector, x, multiplyMatrixVector);
    return vectorResult;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    double startTime, endTime;
    startTime = MPI_Wtime();

    int rank = 0, cntProcess = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    dynamicVector vector = IterativeMethod(rank, cntProcess);

    endTime = MPI_Wtime();

    if(rank == 0) {
        printf("RESULT VECTOR\n");
        PrintVector(vector, SIZE_VECTOR);
        std::cout << "\nTIME: " << endTime - startTime << std::endl;
    }

    delete[] vector;
    MPI_Finalize();
    return 0;
}
