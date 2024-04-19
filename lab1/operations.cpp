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


