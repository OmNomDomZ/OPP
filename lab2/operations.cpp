#include <iostream>
#include <cmath>
#include <omp.h>
#include "operations.h"

void PrintMat(float* mat, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << mat[i * size + j] << "\t";
        }
        std::cout << "\n";
    }
}

void PrintVect(float* vect, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << vect[i] << " ";
    }
}

void FillZero(float*& array, int size) {

    for (int i = 0; i < size; ++i){
        array[i] = 0;
    }

}

void FillMatrix(float*& matrix, int size) {

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[i * size + j] = (i == j) ? 2.0 : 1.0;
        }
    }

}

void FillVectU(float*& u, int size) {

    for (int i = 0; i < size; ++i){
        u[i] = sin((2 * M_PI * i) / size);
    }
}


void FillVectB(float*& b, float*& u, float*& matrix, int size) {
    for (int i = 0; i < size; ++i) {
        float temp = 0.0;
        for (int j = 0; j < size; ++j) {
            temp += matrix[i * size + j] * u[j];
        }
        b[i] += temp;
    }
}
