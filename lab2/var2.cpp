#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <fstream>
#include <thread>
#include <vector>
#include <omp.h>

#include "operations.h"

using std::vector;

const std::size_t VEC_SIZE = 300;
const std::size_t MAT_SIZE = 90000;
const int NUM_THREADS_MAX = omp_get_max_threads();

double* Solve(double*& A, double*& b, int size) {
    double* x = new double[size];
    double* final_x = new double[size];
    FillZero(x, size);
    FillZero(final_x, size);


    double res = EPSILON;
    double norm_x_squared = 0;
    double norm_b = 0;
    for (int i = 0; i < size; ++i)
    {
        norm_b = b[i] * b[i];
    }
    norm_b = sqrt(norm_b);
    int iter = 0;

    #pragma omp parallel
    {
        while (res >= EPSILON) {
            #pragma omp for reduction(+:norm_x_squared)
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    final_x[i] += A[i * size + j] * x[j];           // final_x = Ax
                }
                final_x[i] -= b[i];                                 // final_x = Ax - b
                norm_x_squared += final_x[i] * final_x[i];          // ||Ax - b|| ^ 2
                final_x[i] *= THAU;                                 // final_x = t(Ax - b)
                final_x[i] = x[i] - final_x[i];                     // final_x = x - t(Ax - b)
            }

            #pragma omp for
            for (int i = 0; i < size; i++) {
                x[i] = final_x[i];
            }

            FillZero(final_x, size);

            #pragma omp single
            {
                res = sqrt(norm_x_squared) / norm_b;
                norm_x_squared = 0;
                iter++;
            }

        }
    }

    std::cout << iter << "\n";

    delete[] final_x;
    return x;
}

int main(int argc, char** argv) {

    double* A = new double[MAT_SIZE];
    double* u = new double[VEC_SIZE];
    double* b = new double[VEC_SIZE];

    FillMatrix(A, VEC_SIZE);
    FillVectU(u, VEC_SIZE);
    FillVectB(b, u, A, VEC_SIZE);

    for (int numThreads = 1; numThreads <= NUM_THREADS_MAX; ++numThreads)
    {

        omp_set_num_threads(numThreads);

        std::chrono::high_resolution_clock clock;

        auto start = clock.now();
        double* x = Solve(A, b, VEC_SIZE);
        auto end = clock.now();

        auto time = std::chrono::duration_cast<std::chrono::seconds> (end - start);
        std::cout << time.count() << " s\n";

        delete[] x;
    }

    delete[] A;
    delete[] u;
    delete[] b;
}
