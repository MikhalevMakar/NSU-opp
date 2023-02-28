//
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <assert.h>
//
//enum { SIZE_MATRIX = 1024,
//    ARBITRARY_VALUE = 0,
//    ZERO_VALUE = 0
//};
//
//const double tau =  1e-5;
//const double epsilon = 1e-5;
//
//typedef double* dynamicArray;
//typedef double* dynamicMatrix;
//
//void PrintVector(dynamicArray vector, int size) {
//    for(int j = 0; j < size; ++j) {
//        std::cout << (double)vector[j] << " ";
//    }
//    printf("\n");
//}
//
//dynamicArray GenerateDynamicArray(int size) {
//    dynamicMatrix vector = new double[size];
//    assert(vector != NULL);
//    return vector;
//}
//
//void GenerateVectorArbitraryValue(double* vector, int size, double value) {
//    #pragma omp parallel for
//    for(int i = 0; i < size; ++i) {
//        vector[i] = value;
//    }
//}
//
//dynamicArray GenerateSolutionVector(int newSize) {
//    dynamicArray vector =  GenerateDynamicArray(newSize);
//
//    for(int i = 0; i < newSize; ++i) {
//        vector[i] = (i < SIZE_MATRIX) ? ARBITRARY_VALUE : ZERO_VALUE;
//    }
//    return vector;
//}
//
//dynamicArray GenerateVectorRightParts(const int newSize) {
//    dynamicArray vector =  GenerateDynamicArray(newSize);
//
//    for(int i = 0; i < newSize; ++i) {
//        vector[i] = (i < SIZE_MATRIX) ? SIZE_MATRIX+1 : ZERO_VALUE;
//    }
//    return vector;
//}
//
//int GetBalanceSizeVector(const int& sizeOrigin, const int& countProcess) {
//    int buildNewSize = sizeOrigin;
//    while(buildNewSize % countProcess != 0) {
//        ++buildNewSize;
//    }
//    return buildNewSize;
//}
//
//int GetCntCurrentFillLineMatrix(const int& rank, const int& cntProcess, int newSize) {
//    int resultFillLine = 0;
//    for (int startRank = 1; startRank < rank+1; ++startRank) {
//        resultFillLine += (newSize / cntProcess);
//    }
//    return resultFillLine;
//}
//
//dynamicMatrix GeneratePartMatrix(const int& rank, const int& cntProcess, int fictitiousSize) {
//    int partFictitiousSizeMatrix = fictitiousSize * (fictitiousSize / cntProcess);
//
//    dynamicMatrix partMatrix = GenerateDynamicArray(partFictitiousSizeMatrix);
//
//    int countRows = fictitiousSize / cntProcess;
//    int index = countRows * rank;
//    int  numberCntLine = GetCntCurrentFillLineMatrix(rank, cntProcess,  fictitiousSize);
//
//    for (int i = 0, offset = 0; i < fictitiousSize/cntProcess; ++i, offset += fictitiousSize) {
//        for (int j = 0; j < fictitiousSize; ++j) {
//            partMatrix[offset + j] = (numberCntLine+i+1 <= SIZE_MATRIX) ? 1.0 : ZERO_VALUE;
//        }
//
//        if(numberCntLine+i+1 <= SIZE_MATRIX)
//            partMatrix[index + i] = 2.0;
//
//        index += fictitiousSize;
//    }
//
//    return partMatrix;
//}
//
//dynamicArray MultiplyVectors(const dynamicArray vector1, const dynamicArray vector2,
//                             dynamicArray result, const int& size, const int& cntProcess) {
//    GenerateVectorArbitraryValue(result, size / cntProcess, ZERO_VALUE);
//    for(int i = 0; i < size / cntProcess; ++i) {
//        for  (int j = 0; j < size; ++j) {
//            result[i] += vector1[j+i*size] * vector2[j];
//        }
//    }
//    return result;
//}
//
//dynamicArray MinusVectors(const dynamicArray vector1, const dynamicArray vector2,
//                          dynamicArray result, const int& size) {
//    for(int i = 0; i < size; ++i) {
//        result[i] = vector1[i] - vector2[i];
//    }
//    return result;
//}
//
//dynamicArray MultiplyVectorByConstant(const dynamicArray vector, const double& constant,
//                                      dynamicArray result, const int& size) {
//    for(int i = 0; i < size; ++i) {
//        result[i] = vector[i] * constant;
//    }
//    return result;
//}
//
//double FormingFirstNorm(const dynamicArray vector) {
//    double sumVector = 0;
//    for(int i = 0; i < SIZE_MATRIX; ++i) {
//        sumVector += vector[i]*vector[i];
//    }
//    return sqrt(sumVector);
//}
//
//bool IsFirstNormMoreEpsilon(const double& v1, const double& v2) {
//    return !((v1/v2) < epsilon);
//}
//
//void CopyVector(dynamicArray copyVector, const dynamicArray sourceVector, int size) {
//    for(int i = 0; i < size; ++i) {
//        copyVector [i] = sourceVector[i];
//    }
//}
//
//void  DeleteVectors(dynamicMatrix v1, dynamicArray v2, dynamicArray v3,
//                    dynamicArray v4, dynamicArray v5, dynamicArray v6) {
//    delete[] v1;
//    delete[] v2;
//    delete[] v3;
//    delete[] v4;
//    delete[] v5;
//    delete[] v6;
//}
//
////x^(n+1) = x^n – tau(Ax^n – b)
//double* IterativeMethod(const int& rank, const int& cntProcess) {
//    int fictitiousSize = GetBalanceSizeVector(SIZE_MATRIX, cntProcess);
//    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess, fictitiousSize);
//
//    dynamicArray b = GenerateVectorRightParts(fictitiousSize);
//    dynamicArray x = GenerateSolutionVector(fictitiousSize);
//
//    dynamicMatrix vectorResult = GenerateDynamicArray(fictitiousSize);
//    dynamicMatrix multiplyPartMatrix = GenerateDynamicArray(fictitiousSize / cntProcess);
//    dynamicMatrix vectorUtility = GenerateDynamicArray(fictitiousSize);
//    dynamicMatrix compareVectors = GenerateDynamicArray(fictitiousSize);
//
//    double normB = FormingFirstNorm(b);
//
//    do {
//        multiplyPartMatrix = MultiplyVectors(A, x, multiplyPartMatrix, fictitiousSize, cntProcess);
//
////        MPI_Allgather(multiplyPartMatrix, fictitiousSize / cntProcess, MPI_DOUBLE, vectorUtility,
////                      fictitiousSize / cntProcess, MPI_DOUBLE, MPI_COMM_WORLD);
//
//        compareVectors = MinusVectors(vectorUtility, b, compareVectors, fictitiousSize);
//
//        vectorResult = MinusVectors(x,
//                                    MultiplyVectorByConstant(compareVectors, tau, vectorResult, fictitiousSize),
//                                    vectorResult, fictitiousSize);
//
//        CopyVector(x, vectorResult, fictitiousSize);
//
//    } while(IsFirstNormMoreEpsilon(FormingFirstNorm(compareVectors), normB));
//
//    DeleteVectors(A, b, x, multiplyPartMatrix, vectorUtility, compareVectors);
//    return vectorResult;
//}
//
//int main(int argc, char* argv[]) {
////    MPI_Init(&argc, &argv);
//    double startTime, endTime;
////    startTime = MPI_Wtime();
//
//    int rank = 0, cntProcess = 0;
//
////    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
////    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);
//
//    //dynamicArray vector = IterativeMethod(rank, cntProcess);
//
//    //endTime = MPI_Wtime();
//
////    if(rank == 0) {
////        printf("RESULT VECTOR\n");
////        PrintVector(vector, SIZE_MATRIX);
////        std::cout << "\nTIME: " << endTime - startTime << std::endl;
////    }
////    delete[] vector;
//    //MPI_Finalize();
//    return 0;
//}

#include <math.h>
#define N 10000
float x[N];
int main() { int i;
    float k = 2*3.14159265/N;
#pragma omp parallel for
    for (i=0;i<N;i++) x[i]=sinf(k*i);
    return 0;
}