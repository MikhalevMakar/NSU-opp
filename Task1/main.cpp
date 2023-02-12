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

__attribute__((unused)) void PrintMatrix(double* matrix) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        for(int j = 0; j < SIZE_MATRIX; ++j) {
            std::cout << matrix[i* SIZE_MATRIX + j];
        }
        std::cout << std::endl;
    }
}

void PrintVector(double* vector) {
    for(int j = 0; j < SIZE_MATRIX; ++j) {
        std::cout << (double)vector[j] << " ";
    }
    printf("\n");
}

double* GenerateSolutionVector() {
    double* partMatrix = new double[SIZE_MATRIX];
    assert(partMatrix != NULL);
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        partMatrix[i] = ARBITRARY_VALUE;
    }
    return partMatrix;
}


double* GenerateVectorRightParts() {
    double* vector = new double[SIZE_MATRIX];
    assert(vector != NULL);

    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vector[i] = SIZE_MATRIX+1;
    }
    return vector;
}

double* GenerateVectorArbitrarySizeValue(int size, double value) {
    double* vector = new double[size];
    assert(vector != NULL);

    for(int i = 0; i < size; ++i) {
        vector[i] = value;
    }
    return vector;
}

double* GeneratePartMatrix(const int& rank, const int& countProcess) {
    int partSizeMatrix = SIZE_MATRIX * (SIZE_MATRIX / countProcess);
    double* partMatrix = GenerateVectorArbitrarySizeValue(partSizeMatrix, 1.0);

    for(int  i = 0, offset = 0; i < SIZE_MATRIX / countProcess; ++i) {
        partMatrix[offset+rank+i] = 2.0;
        offset += SIZE_MATRIX;
    }
    return partMatrix;
}

double* MultVectors(const double* vector1, const double* vector2, int cntProcess) {
    double* vectorResult = GenerateVectorArbitrarySizeValue(SIZE_MATRIX, ZERO_VALUE);
    
    for(int i = 0; i < SIZE_MATRIX / cntProcess; ++i) {
        for  (int j = 0; j < SIZE_MATRIX; ++j) {
            vectorResult[i] += vector1[j+i*SIZE_MATRIX] * vector2[j];
        }
    }
    return vectorResult;
}

double* MinusVectors(const double* vector1, const double* vector2) {
    double* vectorResult = new double[SIZE_MATRIX];
    assert(vectorResult != NULL);

    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vectorResult[i] = vector1[i] - vector2[i];
    }
    return vectorResult;

}

double* MultVectorByConstant(double* vector, double constant) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        vector[i] *= constant;
    }
    return vector;
}
double FormingFirstNorm(const double* vector) {
    double sumVector = 0;
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        sumVector += vector[i]*vector[i];
    }
    return sqrt(sumVector);
}

double NormCalculation(const double* multAx,
                       const double* b) {
    return FormingFirstNorm(MinusVectors(multAx, b)) / FormingFirstNorm(b);
}

bool IsFirstNormMoreEpsilon(const double* multAx,
                            const double* b) {
    return !(NormCalculation(multAx, b) < ε);
}

void CopyVector(double* copyVector, const double* sourceVector) {
    for(int i = 0; i < SIZE_MATRIX; ++i) {
        copyVector [i] = sourceVector[i];
    }
}

//x^(n+1) = x^n – τ(Ax^n – b)
void  DeleteVectors(double* v1, double* v2, double* v3) {
    delete[] v1;
    delete[] v2;
    delete[] v3;
}

double* IterativeMethod(int rank, int cntProcess) {
    double* A = GeneratePartMatrix(rank, cntProcess);
    double* b = GenerateVectorRightParts();
    double* x = GenerateSolutionVector();
    //MPI_Barrier(MPI_COMM_WORLD);

     double* vectorResult;
     double* multAx;

    do {
        multAx = MultVectors(A, x, cntProcess);
        double* mult_Ax_minus_b = MinusVectors(multAx, b);
        double* multVectorByConst = MultVectorByConstant(mult_Ax_minus_b, τ);
        vectorResult = MinusVectors(x, multVectorByConst);
        CopyVector(x, vectorResult);

        delete[] multVectorByConst;
        delete[] mult_Ax_minus_b;
    } while(IsFirstNormMoreEpsilon(multAx, b));

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
