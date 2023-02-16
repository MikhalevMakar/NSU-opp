#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>

enum { SIZE_MATRIX = 4,
       ARBITRARY_VALUE = 1,
       ZERO_VALUE = 0
};
//τ
//ε
const double τ =  1e-4;
const double ε = 1e-7;

typedef double* dynamicArray;
typedef double* dynamicMatrix;

void PrintVector(dynamicArray vector, int size) {
    for(int j = 0; j < size; ++j) {
        std::cout << (double)vector[j] << " ";
    }
    printf("\n");
}

void PrintMatrix(dynamicArray vector, int size, int cntProcess) {
    for(int j = 0; j < size/cntProcess; ++j) {
        for(int i = 0; i < size; ++i) {
            std::cout << (double) vector[j*size+i] << " ";
        }
        printf("\n");
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

int GetBalanceSizeVector(const int& sizeOrigin, const int& countProcess) {
    int buildNewSize = sizeOrigin;
    while(buildNewSize % countProcess != 0) {
        ++buildNewSize;
    }
    return buildNewSize;
}

int GetCntFillVectors(const int& rank, const int& cntProcess, int newSize) {
    int resultFillLine = 0;
    for (int startRank = 1; startRank < rank+1; ++startRank) {
        resultFillLine += newSize;
    }
    return resultFillLine;
}


dynamicArray GenerateSolutionVector(int fictitiousSize, const int& rank, const int& cntProcess) {
    dynamicArray vector =  GenerateDynamicArray(fictitiousSize / cntProcess);

    int  numberCntLine = GetCntFillVectors(rank, cntProcess, fictitiousSize / cntProcess);
    int newSize = fictitiousSize/cntProcess;
    for(int i = 0; i < newSize; ++i) {
        vector[i] = (numberCntLine+i < SIZE_MATRIX) ? ARBITRARY_VALUE : ZERO_VALUE;
    }
    return vector;
}

// | 2 1 1 1 0 0 |    1
// | 1 2 1 1 0 0 |    1
// | 1 1 2 1 0 0 |    1
// | 1 1 1 2 0 0 |    1
// | 0 0 0 0 0 0 |    0
// | 0 0 0 0 0 0 |    0

// 5 5 5 5 5
dynamicArray GenerateVectorRightParts(const int& fictitiousSize, const int& rank , const int& cntProcess) {
    dynamicArray vector = GenerateDynamicArray(fictitiousSize/cntProcess);

    int  numberCntLine = GetCntFillVectors(rank, cntProcess, fictitiousSize / cntProcess);
    for(int i = 0; i < fictitiousSize/cntProcess; ++i) {
        vector[i] = (numberCntLine+i < SIZE_MATRIX) ? SIZE_MATRIX + 1 : ZERO_VALUE;
    }
    return vector;
}

int GetCntCurrentFillLineMatrix(const int& rank, const int& cntProcess, int newSize) {
    int resultFillLine = 0;
    for (int startRank = 1; startRank < rank+1; ++startRank) {
        resultFillLine += newSize / cntProcess;
    }
    return resultFillLine;
}

dynamicMatrix GeneratePartMatrix(const int& rank, const int& cntProcess, int fictitiousSize) {
    int partFictitiousSizeMatrix = fictitiousSize * (fictitiousSize / cntProcess);

    dynamicMatrix partMatrix = GenerateDynamicArray(partFictitiousSizeMatrix);

    GenerateVectorArbitraryValue(partMatrix, partFictitiousSizeMatrix, ZERO_VALUE);

    int countRows = fictitiousSize / cntProcess;
    int index = countRows * rank;
    int  numberCntLine = GetCntCurrentFillLineMatrix(rank, cntProcess,  fictitiousSize);

    for (int i = 0, offset = 0; i < fictitiousSize/cntProcess; ++i, offset += fictitiousSize) {
        if(numberCntLine+i+1 <= SIZE_MATRIX) {
            for (int j = 0; j < SIZE_MATRIX; ++j) {
                partMatrix[offset + j] = 1.0;
            }
            partMatrix[index + i] = 2.0;
        }
        index += fictitiousSize;
    }
    return partMatrix;
}

dynamicArray MultiplyVectors(const dynamicArray vector1,
                             const dynamicArray vector2,
                             dynamicArray result,
                             const int& size,
                             const int& rank,
                             const int& cntProcess) {

    GenerateVectorArbitraryValue(result, size/cntProcess, ZERO_VALUE);

    //MPI_Scatter(vector1, size * size / cntProcess, MPI_DOUBLE, );

    dynamicArray  cellArray = GenerateDynamicArray(size/cntProcess);



    for(int i = 0; i < size/cntProcess; ++i) {
        for  (int j = 0; j < size; ++j) {
            //std::cout << vector1[j+i*size] << " " << vector2[i] << "\n";
            result[j] += vector1[j+i*size] * vector2[i];
        }
    }

//    MPI_Allgather(cellArray,
//                  size / cntProcess,
//                  MPI_DOUBLE,
//                  result,
//                  size / cntProcess,
//                  MPI_DOUBLE,
//                  MPI_COMM_WORLD);

    return result;
}

dynamicArray MinusVectors(const dynamicArray vector1, const dynamicArray vector2, dynamicArray result, const int& size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

dynamicArray MultiplyVectorByConstant(const dynamicArray vector, dynamicArray result, const double& constant, const int& size) {
    for(int i = 0; i < size; ++i) {
        result[i] = vector[i] * constant;
    }
    return result;
}

double FormingFirstNorm(const dynamicArray vector) {
    double sumVector = 0;
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        sumVector += vector[i]*vector[i];
    }
    return sqrt(sumVector);
}

double NormCalculation(const dynamicArray vector1, const dynamicArray vector2) {
    return FormingFirstNorm(vector1) / FormingFirstNorm(vector2);
}

bool IsFirstNormMoreEpsilon(const dynamicArray vector1, const dynamicArray vector2) {
    return !(NormCalculation(vector1, vector2) < ε);
}

void CopyVector(dynamicArray copyVector, const dynamicArray sourceVector, int size) {
    for(int i = 0; i < size; ++i) {
        copyVector [i] = sourceVector[i];
    }
}

void  DeleteVectors( dynamicArray v1, dynamicArray v2, dynamicArray v3, dynamicArray v4, dynamicArray v5, dynamicArray v6) {
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
    int fixedCutMatrixSize = fictitiousSize*fictitiousSize/cntProcess;
    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess, fictitiousSize);
    dynamicMatrix final_vect_res = GenerateDynamicArray(fictitiousSize);

    //double* ABuf = new double[fixedCutMatrixSize];


//    MPI_Scatter(A, fixedCutMatrixSize, MPI_DOUBLE, ABuf,
//                fixedCutMatrixSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    dynamicMatrix multiplyPartMatrix = GenerateSolutionVector(fictitiousSize, rank, cntProcess);
    dynamicArray x = GenerateSolutionVector(fictitiousSize, rank, cntProcess);

    //multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, rank, cntProcess);
    //MPI_Allreduce(multiplyPartMatrix, final_vect_res, fictitiousSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //std::cout << "multiply part matrix\n";
    //PrintVector(x, fictitiousSize/cntProcess);

    //std::cout << "Matrix A \n";
    //if(rank == 0) PrintMatrix(A, SIZE_MATRIX, cntProcess);
    //std::cout << "\n_______________________________________\n";

    dynamicArray b = GenerateVectorRightParts(fictitiousSize, rank, cntProcess);

    //std::cout << "vector x\n";
    //PrintVector(x, fictitiousSize/cntProcess);

    //dynamicArray partSolutionVector = GenerateSolutionVector(fictitiousSize, rank, cntProcess);

    //std::cout << "multiply Vector: \n";
    //PrintVector(b, fictitiousSize/cntProcess);
    //std::cout << "\n____________\n";

    dynamicMatrix vectorResult = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix vectorUtility = GenerateDynamicArray(fictitiousSize);
    dynamicMatrix multiplyVectors = GenerateDynamicArray(fictitiousSize);

   // do {

    multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, rank, cntProcess);
//         std::cout << "\nmultVecotrs\n";
        PrintVector(multiplyPartMatrix, fictitiousSize);
//         std::cout << "_________";

//        MPI_Allgather(multiplyPartMatrix,
//                      fictitiousSize / cntProcess,
//                      MPI_DOUBLE,
//                      vectorUtility,
//                      fictitiousSize / cntProcess,
//                      MPI_DOUBLE,
//                      MPI_COMM_WORLD);
//
//        CopyVector(multiplyVectors, vectorUtility, fictitiousSize);
//
//        vectorUtility = MinusVectors(vectorUtility, b, vectorUtility, fictitiousSize);
//
//        vectorResult = MinusVectors(x, MultiplyVectorByConstant(vectorUtility, vectorResult, τ, fictitiousSize), vectorResult, fictitiousSize);
//        CopyVector(x, vectorResult, fictitiousSize);
//    } while(IsFirstNormMoreEpsilon(vectorUtility, b));
//
    DeleteVectors(A, b, vectorUtility, multiplyPartMatrix, multiplyVectors, x);
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

//    if(rank == 0) {
//        printf("RESULT VECTOR\n");
//        PrintVector(vector, SIZE_MATRIX);
//        std::cout << "\nTIME: " << endTime - startTime << std::endl;
//    }

    delete[] vector;
    return 0;
}