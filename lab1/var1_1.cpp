// #include <mpi.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "operations.h"

using namespace std;

const std::size_t VEC_SIZE = 2500;
const std::size_t MAT_SIZE = 6250000;
const float epsilon = 1e-6;
const float tau = 0.001;

using namespace std;

float* readBinaryFile(std::string filename) {
    float* buffer = nullptr;

    std::ifstream file("./data/" + filename, std::ios::binary);

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
    printf("cool\n");
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

void calcPartsMatrix(int* countLines, int* shiftLinesMatrix, int size){
    for (int i = 0; i < size; ++i){
        countLines[i] = VEC_SIZE / size;
    }

    int shift = 0;
    for (int i = 0; i < size; ++i){
        if (i < VEC_SIZE % size){
            countLines[i]++;
        }
        shiftLinesMatrix[i] = shift;
        shift += countLines[i];
    }
}

void calcScattervArgs(int* countLines, int* shiftLinesMatrix, int* countLinesForScatterv, int* shiftLinesMatrixForScatterv, int size){
    for (int i = 0; i < size; ++i){
        countLinesForScatterv[i] = countLines[i] * VEC_SIZE;
        shiftLinesMatrixForScatterv[i] = shiftLinesMatrix[i] * VEC_SIZE;
    }
}


int main(int argc, char **argv) {

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* countLines = new int[size]();
    int* shiftLinesMatrix = new int[size]();

    int* countLinesForScatterv = new int[size]();
    int* shiftLinesMatrixForScatterv = new int[size]();

    calcPartsMatrix(countLines, shiftLinesMatrix, size);

    calcScattervArgs(countLines, shiftLinesMatrix, countLinesForScatterv, shiftLinesMatrixForScatterv, size);

    float* A;
    float* b;
    float* x;

    if (rank == 0){
        A = readBinaryFile("matA.bin");
        b = readBinaryFile("vecB.bin");
        x = new float[VEC_SIZE];
        FillZero(x, VEC_SIZE);
    }

    // рассылаем всем веторы
    MPI_Bcast(b, VEC_SIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, VEC_SIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* PartA = new float[countLines[rank] * VEC_SIZE];
    MPI_Scatterv(A, countLinesForScatterv, shiftLinesMatrixForScatterv, MPI_FLOAT,
                PartA, countLinesForScatterv[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);

    float normVectorAxb = 0.0;
    float normB = 0.0;
    float* AxbVector = new float[countLines[rank]]; // Ax - b
    float* xPart = new float[countLines[rank]];

    // считаем ||b||
    if (rank  == 0){
        for (int i = 0; i < VEC_SIZE; ++i){
            normB += b[i] * b[i];
        }
        normB = sqrt(normB);
    }

    float generalNorm = 0.0;
    float quit = 1;

    double startTime = MPI_Wtime();

    while (quit > epsilon){
        for (int i = 0; i < countLines[rank]; ++i){
            float Ax = 0;
            for (int j = 0; j < VEC_SIZE; ++j){
                Ax += PartA[i * VEC_SIZE + j] * x[j];
            }
            AxbVector[i] = Ax - b[i + shiftLinesMatrix[rank]];
        }

        for (int i = 0; i < countLines[rank]; ++i){
            xPart[i] = x[shiftLinesMatrix[rank] + i] - tau * AxbVector[i]; // x^{n+1} = x^{n} - t*Axb
        }

        for (int  i = 0; i < countLines[rank]; ++i){
            normVectorAxb += AxbVector[i] * AxbVector[i]; // || Ax - b ||
        }

        // собираем части вектора x
        MPI_Allgatherv(xPart, countLines[rank], MPI_FLOAT, x, countLines, shiftLinesMatrix,
                                                                MPI_FLOAT, MPI_COMM_WORLD);

        //собираем все в одно
        MPI_Reduce(&normVectorAxb, &generalNorm, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0){
            quit = sqrt(generalNorm) / normB;
        }
        MPI_Bcast(&quit, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

        normVectorAxb = 0;
        generalNorm = 0;
    }

    double finalTime = MPI_Wtime();
    if (rank == 0){
        cout << "time = " << finalTime - startTime << endl;
    }

    delete[] A;
    delete[] b;
    delete[] x;
    delete[] shiftLinesMatrix;
    delete[] countLines;
    delete[] AxbVector;
    delete[] xPart;

    MPI_Finalize();
    return 0;
}
