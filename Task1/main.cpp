#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

enum { SIZE_MATRIX = 4,
    ARBITRARY_VALUE = 0,
    ZERO_VALUE = 0
};

const double τ =  1e-2;
const double ε = 1e-7;

typedef double* dynamicArray;
typedef double* dynamicMatrix;

__attribute__((unused)) void PrintMatrix(dynamicArray matrix) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        for(int j = 0; j < SIZE_MATRIX; ++j) {
            std::cout << matrix[i* SIZE_MATRIX + j];
        }
        std::cout << std::endl;
    }
}

void PrintVector(dynamicArray vector) {
    for(int j = 0; j < SIZE_MATRIX; ++j) {
        std::cout << (double)vector[j] << " ";
    }
    printf("\n");
}

dynamicArray GenerateSolutionVector() {
    dynamicArray partMatrix = new double[SIZE_MATRIX];
    assert(partMatrix != NULL);
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        partMatrix[i] = ARBITRARY_VALUE;
    }
    return partMatrix;
}


dynamicArray GenerateVectorRightParts() {
    dynamicArray vector = new double[SIZE_MATRIX];
    assert(vector != NULL);

    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vector[i] = SIZE_MATRIX+1;
    }
    return vector;
}

dynamicMatrix GenerateVectorArbitrarySizeValue(int size, double value) {
    dynamicMatrix vector = new double[size];
    assert(vector != NULL);

    for(int i = 0; i < size; ++i) {
        vector[i] = value;
    }
    return vector;
}

dynamicMatrix GeneratePartMatrix(const int& rank, const int& countProcess) {
    int partSizeMatrix = SIZE_MATRIX * (SIZE_MATRIX / countProcess);
    dynamicMatrix partMatrix = GenerateVectorArbitrarySizeValue(partSizeMatrix, 1.0);

    for(int  i = 0, offset = 0; i < SIZE_MATRIX / countProcess; ++i) {
        partMatrix[offset+rank+i] = 2.0;
        offset += SIZE_MATRIX;
    }
    return partMatrix;
}

dynamicArray MultVectors(const double* vector1, const double* vector2, double* result, int cntProcess) {
    dynamicArray vectorResult = GenerateVectorArbitrarySizeValue(SIZE_MATRIX, ZERO_VALUE);
    GenerateVectorArbitrarySizeValue(result, SIZE_MATRIX, ZERO_VALUE);
    for(int i = 0; i < SIZE_MATRIX / cntProcess; ++i) {
        for  (int j = 0; j < SIZE_MATRIX; ++j) {
            vectorResult[i] += vector1[j+i*SIZE_MATRIX] * vector2[j];
        }
    }
    return vectorResult;
}

dynamicArray MinusVectors(const double* vector1, const double* vector2) {
    dynamicArray vectorResult = new double[SIZE_MATRIX];
    assert(vectorResult != NULL);

    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vectorResult[i] = vector1[i] - vector2[i];
    }
    return vectorResult;

}

dynamicArray MultVectorByConstant(dynamicArray vector, double constant) {
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

double NormCalculation(const dynamicArray multAx,
                       const dynamicArray b) {
    return FormingFirstNorm(MinusVectors(multAx, b)) / FormingFirstNorm(b);
}

bool IsFirstNormMoreEpsilon(const dynamicArray multAx,
                            const dynamicArray b) {
    return !(NormCalculation(multAx, b) < ε);
}

void CopyVector(dynamicArray copyVector, const dynamicArray sourceVector) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        copyVector [i] = sourceVector[i];
    }
}

//x^(n+1) = x^n – τ(Ax^n – b)
void  DeleteVectors(dynamicArray v1, dynamicArray v2, dynamicArray* v3) {
    delete[] v1;
    delete[] v2;
    delete[] v3;
}

double* IterativeMethod(int rank, int cntProcess) {
    dynamicMatrix A = GeneratePartMatrix(rank, cntProcess);
    dynamicArray b = GenerateVectorRightParts();
    dynamicArray x = GenerateSolutionVector();
    //MPI_Barrier(MPI_COMM_WORLD);

     double* vectorResult;
     double* multiplyVectors;
     double* vectorUtility;
    do {
        multiplyVectors = MultVectors(A, x, cntProcess);
        double* mult_Ax_minus_b = MinusVectors(multiplyVectors, b);
        double* multVectorByConst = MultVectorByConstant(mult_Ax_minus_b, τ);
        vectorResult = MinusVectors(x, multVectorByConst);
        CopyVector(x, vectorResult);
    } while(IsFirstNormMoreEpsilon(multiplyVectors, b));

    DeleteVectors(A, b, x);
    return vectorResult;
}


int main(int argc, char* argv[]) {
    //if(argc != 1) return 1;

    MPI_Init(&argc, &argv);
    int rank = 0, cntProcess = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &cntProcess);

    double* vector = IterativeMethod(rank, cntProcess);

    //if(rank == cntProcess-1) {
    PrintVector(vector);
    //}
    MPI_Finalize();

    return 0;
}
