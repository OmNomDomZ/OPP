#include <mpi.h>
#include <cstdlib>
#include <string>

#define NUM_OF_DIMS_OF_CARTESIAN_GRID 2

MPI_Comm GridComm; // коммуникатор для сетки
MPI_Comm ColComm; // коммуникатор для столбцов
MPI_Comm RowComm; // коммуникатор для строк

int GridCoords[2];
int ProcNum = 0;
int ProcRank = 0;

const int p1 = 4;
const int p2 = 4;

// количество строк матрицы А
const int n1 = 2000;
// количество столбцов матрицы А / строк матрицы B
const int n2 = 2000;
// количество столбцов матрицы В
const int n3 = 2000;

double randDouble() {
    return (double) rand() / RAND_MAX * 50.0 - 5.0;
}

void dataInit(double *pMatrix, int rowCount, int colCount) {
    for (int i = 0; i < rowCount; i++) {
        for (int j = 0; j < colCount; j++) {
            pMatrix[i * colCount + j] = randDouble();
        }
    }
}

void setToZero(double *pMatrix, int rowCount, int colCount) {
    for (int i = 0; i < rowCount; i++) {
        for (int j = 0; j < colCount; j++) {
            pMatrix[i * colCount + j] = 0;
        }
    }
}

void printVector(const double *pVector, int size, int procNum) {
    printf("proc #%d ", procNum);
    for (int i = 0; i < size; i++)
        printf("%7.4f ", pVector[i]);
    printf("\n");
}

void printVector(const int *pVector, int size, int procNum) {
    printf("proc # %d ", procNum);
    for (int i = 0; i < size; i++)
        printf("%d ", pVector[i]);
    printf("\n");
}

void printMatrix(const double *pMatrix, int rowCount, int colCount) {
    for (int i = 0; i < rowCount; i++) {
        for (int j = 0; j < colCount; j++)
            printf("%7.4f ", pMatrix[i * colCount + j]);
        printf("\n");
    }
}

void matrixMul(const double *pAMatrix, const double *pBMatrix, double *pCMatrix, int heightA, int widthAheightB,
               int widthB) {
    for (int i = 0; i < heightA; i++) {
        for (int j = 0; j < widthB; j++)
            for (int k = 0; k < widthAheightB; k++)
                pCMatrix[i * widthB + j] += pAMatrix[i * widthAheightB + k] * pBMatrix[k * widthB + j];
    }
}

void createGridCommunicators() {
    // Количество процессов в каждом измерении сетки
    int dimSize[NUM_OF_DIMS_OF_CARTESIAN_GRID];

    // сигналим о том, переодические у нас или нет
    // Периодическая означает, что края сетки соединены
    int periodic[NUM_OF_DIMS_OF_CARTESIAN_GRID];

    // фиксированный ли у нас размер сетки или нет
    int subDimension[NUM_OF_DIMS_OF_CARTESIAN_GRID];

    dimSize[0] = p1;
    dimSize[1] = p2;

    periodic[0] = 0;
    periodic[1] = 0;

    // создаем коммуникатор
    // reorder 1 - может переупорядочить процессы, оптимизирует
    MPI_Cart_create(MPI_COMM_WORLD, NUM_OF_DIMS_OF_CARTESIAN_GRID, dimSize, periodic, 1, &GridComm);

    // определеняем координаты для каждого процесса
    // принимает номер процесса и возвращает его координаты в решетке
    MPI_Cart_coords(GridComm, ProcRank, NUM_OF_DIMS_OF_CARTESIAN_GRID, GridCoords);

    // нам нужно x и y, поэтому сохраняем лишь одно измерение. 1 - сохр, 0 - несохр
    subDimension[0] = 0;
    subDimension[1] = 1;

    // создает новый коммуникатор, который содержит только подгруппу процессов исходного коммуникатора
    // каждая строка сетки отвечает за n1 / p1 строк матрицы
    MPI_Cart_sub(GridComm, subDimension, &RowComm);

    subDimension[0] = 1;
    subDimension[1] = 0;

    // создаем новую подгруппу колон
    MPI_Cart_sub(GridComm, subDimension, &ColComm);
}

void init(double *&pAMatrix, double *&pBMatrix, double *&pCMatrix,
          double *&pAblock, double *&pBblock, double *&pCblock, int &ABlockSize, int &BBlockSize) {

    ABlockSize = n1 / p1;
    BBlockSize = n3 / p2;

    pAblock = new double[n2 * ABlockSize];
    pBblock = new double[n2 * BBlockSize];
    pCblock = new double[ABlockSize * BBlockSize];

    if (ProcRank == 0) {
        pAMatrix = new double[n1 * n2];
        pBMatrix = new double[n2 * n3];
        pCMatrix = new double[n1 * n3];
        dataInit(pAMatrix, n1, n2);
        dataInit(pBMatrix, n2, n3);
        setToZero(pCMatrix, n1, n3);
    }

    setToZero(pCblock, ABlockSize, BBlockSize);
}

void terminate(const double *AMatrix, const double *BMatrix,
               const double *CMatrix, const double *Ablock, const double *Bblock, const double *Cblock) {

    if (ProcRank == 0) {
        delete[] AMatrix;
        delete[] BMatrix;
        delete[] CMatrix;
    }
    delete[] Ablock;
    delete[] Bblock;
    delete[] Cblock;
}

void dataDistribution(double *AMatrix, double *BMatrix, double *Ablock,
                      double *Bblock, int ABlockSize, int BBlockSize) {

    // только процесс с ранго   м 0 в каждой колонке будет отправлять данные
    // MPI_Scatter разбрасывает данные из массива AMatrix по всем процессам, принадлежащим коммуникатору ColComm
    if (GridCoords[1] == 0) {
        // разрезаем по OY
        MPI_Scatter(AMatrix, ABlockSize * n2, MPI_DOUBLE, Ablock,
                    ABlockSize * n2, MPI_DOUBLE, 0, ColComm);
    }

    // лютейшая рассылочка вдоль OX
    // каждый процесс получает матрицу ABlock [n1 / p1 * n2]
    // MPI_Bcast транслирует данные из блока Ablock процесса с рангом 0 всем остальным процессам в коммуникаторе RowComm
    MPI_Bcast(Ablock, ABlockSize * n2, MPI_DOUBLE, 0, RowComm);


    // распределения данных для матрицы B
    // производные типы
    MPI_Datatype col, colType;

    // создаем дискриптор блока для B
    // тип данных представляет собой столбец матрицы B, с параметрами:
    // n2: количество элементов в столбце,
    // BBlockSize: количество элементов в блоке,
    // n3: смещение между столбцами в матрице.
    MPI_Type_vector(n2, BBlockSize, n3, MPI_DOUBLE, &col);
    // зафиксировали
    MPI_Type_commit(&col);

    // создаем единую колонку из блока
    // BBlockSize * sizeof(double) -- новый размер блока данных
    // передаем сразу блок матрицы В
    MPI_Type_create_resized(col, 0, BBlockSize * sizeof(double), &colType);
    MPI_Type_commit(&colType);

    if (GridCoords[0] == 0) {
        // разрезаем по OX
        MPI_Scatter(BMatrix, 1, colType, Bblock, n2 * BBlockSize,
                    MPI_DOUBLE, 0, RowComm);
    }

    // рассылочка вдоль OY
    MPI_Bcast(Bblock, BBlockSize * n2, MPI_DOUBLE, 0, ColComm);
}

int main(int argc, char *argv[]) {
    double *AMatrix;
    double *BMatrix;
    double *CMatrix;

    int ABlockSize = 0;
    int BBlockSize = 0;

    // текущие блоки матриц
    double *Ablock;
    double *Bblock;
    double *Cblock;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    createGridCommunicators();

    init(AMatrix, BMatrix, CMatrix, Ablock, Bblock,
         Cblock, ABlockSize, BBlockSize);

    double startTime;

    if (ProcRank == 0) {
        startTime = MPI_Wtime();
//        printf("Matrix A \n");
//        printMatrix(AMatrix, n1, n2);
//        printf("Matrix B \n");
//        printMatrix(BMatrix, n2, n3);
    }

    // распределяем данные по блокам
    dataDistribution(AMatrix, BMatrix, Ablock, Bblock, ABlockSize, BBlockSize);

    matrixMul(Ablock, Bblock, Cblock, ABlockSize, n2, BBlockSize);

    // для сбора создаем новый
    MPI_Datatype block, blockType;
    MPI_Type_vector(ABlockSize, BBlockSize, n3, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);

    MPI_Type_create_resized(block, 0, BBlockSize * sizeof(double), &blockType);
    MPI_Type_commit(&blockType);

    // считаем смещение
    int *displ = new int[p1 * p2];
    // количество элементов от каждого процесса
    // каждый элемент в массиве соответствует рангу процесса отправки.
    int *rcount = new int[p1 * p2];
    int blockCount = 0;
    int blockSize = ABlockSize * BBlockSize;
    int numCount = 0;
    int j = 0;

    while (numCount < p1 * p2 * blockSize) {
        int written = 0;
        for (int i = 0; i < n3; i += BBlockSize) {
            displ[j] = blockCount;
            rcount[j] = 1;
            j++;
            blockCount++;
    
            written++;
        }
        numCount += written * blockSize;
        blockCount += written * (ABlockSize - 1);
    }

    MPI_Gatherv(Cblock, blockSize, MPI_DOUBLE, CMatrix,
                rcount, displ, blockType, 0, MPI_COMM_WORLD);

    if (ProcRank == 0) {
        double endTime = MPI_Wtime();
        printf("Matrix C \n");

        printMatrix(CMatrix, n1, n3);
        printf("That took %lf seconds\n", endTime - startTime);
    }

    terminate(AMatrix, BMatrix, CMatrix, Ablock, Bblock, Cblock);
    delete[] displ;
    delete[] rcount;
    MPI_Finalize();

    return 0;
}
