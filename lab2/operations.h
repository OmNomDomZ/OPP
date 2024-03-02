#ifndef OPERATIONS_H_INLUDED
#define OPERATIONS_H_INLUDED

#define THAU 0.0001
#define EPSILON 1.0e-6

void PrintMat(double* mat, int size);

void PrintVect(double* vect, int size);

void FillZero(double*& array, int size);

void FillMatrix(double*& matrix, int size);

void FillVectU(double*& u, int size);

void FillVectB(double*& b, double*& u, double*& matrix, int size);

#endif
