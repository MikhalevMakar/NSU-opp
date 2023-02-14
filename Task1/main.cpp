#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

enum { SIZE_MATRIX = 8,
    ARBITRARY_VALUE = 1,
    ZERO_VALUE = 0
};

const double τ =  1e-2;
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

void GenerateVectorArbitraryValue(double* vector, int size, double value) {
    for(int i = 0; i < size; ++i) {
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
// 2 1 1 1 0 0
// 1 2 1 1 0 0

// 1 1 2 1 0 0
// 1 1 1 2 0 0

// 0 0 0 0 0 0
// 0 0 0 0 0 0

dynamicMatrix GeneratePartMatrix(const int& rank, const int& countProcess, int fictitiousSize) {
    int partFictitiousSizeMatrix = fictitiousSize * (fictitiousSize / countProcess);
    int partOriginSizeMatrix = SIZE_MATRIX * (SIZE_MATRIX / countProcess);

    dynamicMatrix partMatrix = GenerateDynamicArray(partFictitiousSizeMatrix);

    GenerateVectorArbitraryValue(partMatrix, partFictitiousSizeMatrix, ZERO_VALUE);

    int countRows = fictitiousSize / countProcess;
    int startPosition = countRows * rank;
    for (int i = 0, offset = 0; i < fictitiousSize/countProcess && rank < SIZE_MATRIX/countProcess; ++i,  offset += fictitiousSize) {
        for (int j = 0; j < SIZE_MATRIX; ++j) {
            partMatrix[j+offset] = 1.0;
        }
        partMatrix[startPosition + i] = 2.0;
        startPosition += fictitiousSize;
    }
    return partMatrix;
}

//2 1 1 1 0 0 1 2 1 1 0 0

dynamicArray MultiplyVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, int size, int cntProcess) {
    GenerateVectorArbitraryValue(result, size/cntProcess, ZERO_VALUE);
    for(int i = 0; i < size / cntProcess; ++i) {
        for  (int j = 0; j < size; ++j) {
            result[i] += vector1[j+i*size] * vector2[j];
        }
    }
    return result;
}

dynamicArray MinusVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, int size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

dynamicArray MultiplyVectorByConstant(dynamicArray vector, double constant, int size) {
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

double NormCalculation(const dynamicArray vector1, const dynamicArray vector2, dynamicArray vectorUtility, int size) {
    return FormingFirstNorm(MinusVectors(vector1, vector2, vectorUtility, size)) / FormingFirstNorm(vector2);
}

bool IsFirstNormMoreEpsilon(const dynamicArray vector1, const dynamicArray vector2, dynamicArray vectorUtility, int size) {
    return !(NormCalculation(vector1, vector2, vectorUtility, size) < ε);
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


// 2 1 1 1 0 0
// 1 2 1 1 0 0

// 1 1 2 1 0 0
// 1 1 1 2 0 0

// 0 0 0 0 0 0
// 0 0 0 0 0 0

//x^(n+1) = x^n – τ(Ax^n – b)
double* IterativeMethod(const int& rank, const int& cntProcess) {
    int fictitiousSize = GetBalanceSizeVector(SIZE_MATRIX, cntProcess);
    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess, fictitiousSize);

    std::cout<< fictitiousSize << "\n";
    PrintVector(A, fictitiousSize*fictitiousSize/ cntProcess);
    std::cout << "\n";

    dynamicArray b = GenerateVectorRightParts(fictitiousSize, SIZE_MATRIX);
    dynamicArray x = GenerateSolutionVector(fictitiousSize, SIZE_MATRIX);

//    PrintVector(b, fictitiousSize);
//    PrintVector(x, fictitiousSize);

    dynamicMatrix vectorResult = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix multiplyPartMatrix = GenerateDynamicArray(fictitiousSize / cntProcess);
    dynamicMatrix vectorUtility = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix multiplyVectors = GenerateDynamicArray(fictitiousSize);

    int cnt = 0;
    do {
        //PrintVector(A, fictitiousSize / cntProcess);
          multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, cntProcess);
//        std::cout << "rank: " << rank << std::endl;
//        PrintVector(multiplyPartMatrix, 2);
//        std::cout  << std::endl;

//        printf("multiplyPartMatrix: \n ");
//        PrintVector(multiplyPartMatrix, fictitiousSize/cntProcess);
//        printf("\n");
        MPI_Allgather(multiplyPartMatrix,
                      fictitiousSize / cntProcess,
                      MPI_DOUBLE,
                      vectorUtility,
                      fictitiousSize / cntProcess,
                      MPI_DOUBLE,
                      MPI_COMM_WORLD);

        CopyVector(multiplyVectors, vectorUtility, fictitiousSize);
//        if(rank == 1 && cnt==1) PrintVector(vectorResult, 6);

        vectorResult = MinusVectors(vectorUtility, b, vectorResult, fictitiousSize);
//        if(rank == 1 && cnt==1) {
//            PrintVector(vectorUtility, 6);
//            //PrintVector(b, 6);
//            //PrintVector(vectorResult, 6);
//        }

        vectorResult = MultiplyVectorByConstant(vectorResult, τ, fictitiousSize);
//        //if(rank == 1 && cnt==1) PrintVector(vectorResult, 6);
//        //PrintVector(vectorResult, fictitiousSize);
//
        vectorResult = MinusVectors(x, vectorResult, vectorResult, fictitiousSize);
//        //if(rank == 1 && cnt==1) PrintVector(vectorResult, 6);
        CopyVector(x, vectorResult, fictitiousSize);

        //if(rank == 1 && cnt==1) PrintVector(x, 6);
        } while(IsFirstNormMoreEpsilon(multiplyVectors, b, vectorUtility, fictitiousSize));
        //cnt++;
    //} while(cnt != 3);

    DeleteVectors(A, b, x, multiplyPartMatrix, vectorUtility, multiplyVectors);
    return vectorResult;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank = 0, cntProcess = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    dynamicArray vector = IterativeMethod(rank, cntProcess);

    MPI_Finalize();
    //printf("result\n");
    if(rank == 0) PrintVector(vector, SIZE_MATRIX);
    delete[] vector;
    return 0;
}