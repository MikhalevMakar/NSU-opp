#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <omp.h>

enum {
    SIZE_VECTOR = 25,
    ARBITRARY_VALUE = 0
};

const double tau = 1e-5;
const double epsilon = 1e-6;

typedef double *dynamicVector;
typedef double *dynamicMatrix;

void PrintVector(dynamicVector array, const int size) {
    for (int j = 0; j < size; ++j) {
        std::cout << (double) array[j] << " ";
    }
    printf("\n");
}

dynamicVector GenerateDynamicVector(const int size) {
    auto vector = new double[size];
    assert(vector != NULL);
    return vector;
}

void GenerateVectorArbitraryValue(dynamicVector vector, const int size, const double value) {
    for (int i = 0; i < size; ++i) {
        vector[i] = value;
    }
}

dynamicVector GenerateSolutionVector(const int size) {
    dynamicVector vector = GenerateDynamicVector(size);
    GenerateVectorArbitraryValue(vector, size, ARBITRARY_VALUE);
    return vector;
}

dynamicVector GenerateVectorRightParts(const int size) {
    dynamicVector vector = GenerateDynamicVector(size);
    GenerateVectorArbitraryValue(vector, size, size + 1);
    return vector;
}

dynamicMatrix GenerateMatrix(const int size) {
    dynamicMatrix matrix = GenerateDynamicVector(size * size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[size * i + j] = 1.0f;
        }
        matrix[size * i + i] = 2.0f;
    }
    return matrix;
}

void MultiplyMatrixByVector(const dynamicMatrix matrix, const  dynamicVector vector,
                            dynamicVector result, const int cntLineMatrix) {
    for (int i = 0; i < cntLineMatrix; ++i) {
        for (int j = 0; j < SIZE_VECTOR; ++j) {
            result[i] += matrix[i * SIZE_VECTOR + j] * vector[j];
        }
    }
}

void MinusVectors(const  dynamicVector vector1, const dynamicVector vector2,
                  dynamicVector result, const int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = vector1[i] - vector2[i];
    }
}

void MultiplyVectorByConstant(const  dynamicVector vector, const double constant,
                              dynamicVector result, const int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = vector[i] * constant;
    }
}

double FormingEuclideanNorm(const  dynamicVector vector, const int size) {
    double sumVector = 0.0f;
    for (int i = 0; i < size; ++i) {
        sumVector += vector[i] * vector[i];
    }
    return sumVector;
}

bool IsFirstNormMoreEpsilon(const double v1, const double v2) {
    return !((v1 / v2) < epsilon);
}

void CopyVector(dynamicVector copyVector, const  dynamicVector sourceVector, const int size) {
    for (int i = 0; i < size; ++i) {
        copyVector[i] = sourceVector[i];
    }
}


void DeleteVectors(const dynamicMatrix v1, const dynamicVector v2, const dynamicVector v3,
                   const dynamicVector v4, const dynamicVector v5, const dynamicVector v6) {
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    delete[] v5;
    delete[] v6;
}

std::vector<int> GetCountLine(int countThread, int sizeVector) {
    std::vector<int> vectorCountLine;
    vectorCountLine.resize(sizeVector * sizeVector, sizeVector / countThread);
    int countPartLine = sizeVector % countThread;
    for (int i = 0; i < countThread; ++i) {
        if (i < countPartLine) ++vectorCountLine[i];
    }
    return vectorCountLine;
}

dynamicVector IterativeMethod(const int size) {
    dynamicMatrix A = GenerateMatrix(size);

    dynamicVector b = GenerateVectorRightParts(size);
    dynamicVector x = GenerateSolutionVector(size);

    dynamicMatrix vectorResult = GenerateDynamicVector(size);
    dynamicMatrix multiplyPartMatrix = GenerateDynamicVector(size);
    dynamicMatrix vectorUtility = GenerateDynamicVector(size);
    dynamicMatrix compareVectors = GenerateDynamicVector(size);

    double normB = FormingEuclideanNorm(b, size);
    std::vector<int> vectorCountLine;
    bool run;
    double firstNorm, partFirstNorm;

#pragma omp parallel private(partFirstNorm)
{
    int currentThread = omp_get_thread_num();
    int numThread = omp_get_num_threads();

#pragma omp single
{
    vectorCountLine = GetCountLine(numThread, SIZE_VECTOR);
    for(auto value : vectorCountLine) {
        std::cout << value << " ";
    }

    std::cout << "\n size:  " <<  vectorCountLine.size() << std::endl;
}

    int sizePartVector = size / numThread;
    int offsetMatrix = vectorCountLine[currentThread] * size;
    int offsetVector = currentThread * sizePartVector;

    //do {

            MultiplyMatrixByVector(A + offsetMatrix,
                                   x,
                                   multiplyPartMatrix + offsetVector,
                                   sizePartVector);

#pragma omp barrier
        //if(currentThread == 12) PrintVector(multiplyPartMatrix, size);

        MinusVectors(multiplyPartMatrix + offsetVector,
                     b + offsetVector, compareVectors + offsetVector, sizePartVector);

        MultiplyVectorByConstant(compareVectors + offsetVector, tau, vectorResult + offsetVector, sizePartVector);

        MinusVectors(x + offsetVector,
                     vectorResult + offsetVector,
                     vectorResult + offsetVector,
                     sizePartVector);

        CopyVector(x + offsetVector, vectorResult + offsetVector, sizePartVector);

        partFirstNorm = FormingEuclideanNorm(compareVectors + offsetVector, sizePartVector);

#pragma omp single
        firstNorm = 0.0f;
#pragma omp atomic
        firstNorm += partFirstNorm;

#pragma omp barrier
#pragma omp single
        run = IsFirstNormMoreEpsilon(sqrt(firstNorm), normB);
    //} while (run);
}
    DeleteVectors(A, b, x, multiplyPartMatrix, vectorUtility, compareVectors);
    return vectorResult;
}

int main(int argc, char *argv[]) {

    double startTime = omp_get_wtime();

    dynamicVector vector = IterativeMethod(SIZE_VECTOR);

    double endTime = omp_get_wtime();

    printf("RESULT VECTOR\n");
    PrintVector(vector, SIZE_VECTOR);
    std::cout << std::endl << "TIME: " << endTime - startTime << std::endl;

    delete[] vector;
    return 0;
}