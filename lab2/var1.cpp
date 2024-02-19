#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Nx 3
#define Ny 3
#define EPSILON 1e-5  // Погрешность

double* SimpleIterationMethod(double** A, double* prevX, double* b, int N) {
    double t = 0.01;
    double* nextX = (double*)malloc(N * sizeof(double));
    double* mult = (double*)calloc(N, sizeof(double));


    double numerator = 0.0;
    double denominator = 0.0;

    while ((numerator/denominator) < EPSILON){

        for (int i = 0; i < N; ++i){
            for (int j = 0; j < N; ++j){
                mult[i] += A[i][j] * prevX[j];
            }
        }

        for (int i = 0; i < N; ++i){
            numerator += pow(mult[i] - b[i], 2);
        }
        numerator = sqrt(numerator);

        for (int i = 0; i < N; ++i){
            denominator += pow(b[i], 2);
        }
        denominator = sqrt(denominator);


        for (int i = 0; i < N; ++i){
            nextX[i] = prevX[i] - t * (mult[i] - b[i]);
        }

        for (int i = 0; i < N; ++i){
            prevX[i] = nextX[i];
        }
    }

    free(mult);
    return nextX;
}

int main() {

    int N = 3;

    // double* x0 = (double*)malloc (N * sizeof(double));
    double x0[] = {1, 2, 3};
    double b[] = {1, 2, 3};
    // double* b = (double*)malloc (N * sizeof(double));
    double a = 0.0;

    double** A = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; ++i) {
        A[i] = (double*)malloc(N * sizeof(double));
        for (int j = 0; j < N; ++j) {
            A[i][j] = a++;
        }
    }

    double* aaa = SimpleIterationMethod(A, x0, b, N);

    for (int i = 0; i < N; ++i)
    {
        printf("%f\n", aaa[i]);
    }

    for (int i = 0; i < N; ++i) {
        free(A[i]);
    }
    free(A);

    return 0;
}
