#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <omp.h>

#include "operations.h"

const std::size_t VEC_SIZE = 2500;
const std::size_t MAT_SIZE = 6250000;
const int NUM_THREADS_MAX = omp_get_max_threads();

float* readBinaryFile(const char* filename) {
    float* buffer = nullptr;

    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return nullptr;
    }

    file.seekg(0, std::ios::end);
    std::size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    buffer = new float[fileSize / sizeof(float)];

    file.read(reinterpret_cast<char*>(buffer), fileSize);

    file.close();
    return buffer;
}

void writeFloatArrayToFile(const char* filename, float* array, size_t size) {
    FILE* file = fopen(filename, "wb");

    if (file != NULL) {
        fwrite(array, sizeof(float), size, file);

        fclose(file);
    } else {
        printf("Не удалось открыть файл для записи\n");
    }
}

float* Solve(float*& A, float*& b, float* x, int size) {
    float* final_x = new float[size];
    FillZero(final_x, size);

    float res = 1;
    float norm_x_squared = 0;
    float norm_b = 0;
    for (int i = 0; i < size; ++i)
    {
        norm_b += b[i] * b[i];
    }
    norm_b = sqrt(norm_b);

    int iter = 0;

    while (res >= EPSILON) {
        #pragma omp parallel for reduction(+:norm_x_squared)
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                final_x[i] += A[i * size + j] * x[j];
            }
            final_x[i] -= b[i];
            norm_x_squared += final_x[i] * final_x[i];
            final_x[i] = final_x[i] * THAU;
            final_x[i] = x[i] - final_x[i];
        }

        #pragma omp parallel for
        for (int i = 0; i < size; i++) {
            x[i] = final_x[i];
        }

        res = sqrt(norm_x_squared) / norm_b;
        norm_x_squared = 0;
        FillZero(final_x, size);
        iter++;
    }
    std::cout << iter << "\n";

    delete[] final_x;
    return x;
}

int main(int argc, char** argv) {

    float* A;
    float* b;
    float* x;

    A = readBinaryFile("matA.bin");
    b = readBinaryFile("vecB.bin");
    x = new float[VEC_SIZE];

    for (int numThreads = 1; numThreads <= NUM_THREADS_MAX; ++numThreads)
    {
        FillZero(x, VEC_SIZE);
        omp_set_num_threads(numThreads);

        std::chrono::high_resolution_clock clock;

        auto start = clock.now();

        x = Solve(A, b, x, VEC_SIZE);

        auto end = clock.now();

        auto time = std::chrono::duration_cast<std::chrono::seconds> (end - start);
        std::cout << time.count() << " s\n";
    }

    writeFloatArrayToFile("myVecX", x, VEC_SIZE);

    delete[] x;
    delete[] A;
    delete[] b;
}
