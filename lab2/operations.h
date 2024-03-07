#ifndef OPERATIONS_H_INLUDED
#define OPERATIONS_H_INLUDED

#define THAU -0.01
#define EPSILON 1.0e-3

void PrintMat(float* mat, int size);

void PrintVect(float* vect, int size);

void FillZero(float*& array, int size);

void FillMatrix(float*& matrix, int size);

void FillVectU(float*& u, int size);

void FillVectB(float*& b, float*& u, float*& matrix, int size);

#endif
