#include <iostream>
#include <cstdlib>
#include "mpi.h"
using namespace std;

int main(int argc, char* argv[])
{
    int ProcNum, ProcRank;
    int Rows, Cols;
    MPI_Status status;
    if (argc < 3)
        return 1;
    Rows = atoi(argv[2]);
    Cols = atoi(argv[1]);
    if ( Cols < 1 || Rows < 1)
        return 1;
    if (argc - 3 != Rows*Cols)
        return 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcNum == 1)
    {
        double** array = new double*[Rows];
        for (int i = 0; i < Rows; i++)
            array[i] = new double[Cols];
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Cols; j++)
                array[i][j] = atof(argv[j + Cols*i + 3]);
        for (int i = 0; i < Rows; i++)
        {
            int max = 0;
            for (int j = 1; j < Cols; j++)
            {
                if (array[i][j] > array[i][max])
                    max = j;
            }
            cout << i + 1 << " Row's maximum is " << max + 1 << " element with " << array[i][max] << " value" << endl;
        }
        for (int i = 0; i < Rows; i++)
            delete[] array[i];
        delete[] array;
    }
    else
    {
        if (ProcRank == 0)
        {
            double** array = new double*[Rows];
            for (int i = 0; i < Rows; i++)
                array[i] = new double[Cols];
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Cols; j++)
                    array[i][j] = atof(argv[j + Cols*i + 3]);

            int RowPerProc = Rows / ProcNum;
            int MainProcRows = RowPerProc + (Rows % ProcNum);
            MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 1; i < ProcNum; i++)
                for (int j = 0; j < RowPerProc; j++)
                {
                    MPI_Send(&(array[(i - 1)*RowPerProc + j][0]), Cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                }
            for (int i = (ProcNum-1)*RowPerProc; i < Rows; i++)
            {
                int max = 0;
                for (int j = 1; j < Cols; j++)
                {
                    if (array[i][j] > array[i][max])
                        max = j;
                }
                cout << i + 1 << " Row's maximum is " << max + 1 << " element with " << array[i][max] << " value" << endl;
            }
            for (int i = 0; i < Rows; i++)
                delete[] array[i];
            delete[] array;
        }
        else
        {
            int RowPerProc;
            MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
            //        MPI_Recv(&RowPerProc, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (RowPerProc > 0)
            {
                double** subarray = new double*[RowPerProc];
                for (int i = 0; i < RowPerProc; i++)
                    subarray[i] = new double[Cols];
                for (int i = 0; i < RowPerProc; i++)
                {
                    MPI_Recv(&subarray[i][0], Cols, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }
                for (int i = 0; i < RowPerProc; i++)
                {
                    int max = 0;
                    for (int j = 0; j < Cols; j++)
                    {
                        if (subarray[i][j] > subarray[i][max])
                            max = j;
                    }
                    cout << i + (ProcRank - 1)*RowPerProc + 1 << " Row's maximum is " << max + 1 << " element with " << subarray[i][max] << " value" << endl;
                }
                for (int i = 0; i < RowPerProc; i++)
                    delete[] subarray[i];
                delete[] subarray;
            }
        }

    }

//    cout << "Hello, world";
    MPI_Finalize();

    return 0;
}
