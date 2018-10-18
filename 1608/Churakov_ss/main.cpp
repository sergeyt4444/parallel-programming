#include <iostream>
#include "time.h"
#include <cstdlib>
#include "mpi.h"
using namespace std;

void RandomizeArray(double* rarray, int size)
{
	for (int i = 0; i < size; i++)
	{
		*(rarray + i) = rand() % 101 - 50;
	}
}

void MatrixFindMaxPerRow(double* matrix, int Rows, int Cols, int* maxarray)
{
	for (int i = 0; i < Rows; i++)          
	{
		maxarray[i] = 0;
		double maxElem = matrix[i*Cols];
		for (int j = 1; j < Cols; j++)
		{
			if (matrix[i*Cols + j] > maxElem)
			{
				maxElem = matrix[i*Cols + j];
				maxarray[i] = j;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	srand(time(0));
    int ProcNum, ProcRank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcNum == 1)
    {
		double timef, times, tickt;
		timef = MPI_Wtime();
		tickt = MPI_Wtick();
		int Rows, Cols; 
		if (argc != 3)
			return 1;
		Rows = atoi(argv[2]);  
		Cols = atoi(argv[1]);
		if (Cols < 1 || Rows < 1)
			return 1;
		double* array = new double[Rows*Cols];   
		RandomizeArray(array, Rows*Cols);
		int* maxarray = new int[Rows];
		MatrixFindMaxPerRow(array, Rows, Cols, maxarray);
		for (int i = 0; i < Rows; i++)
			cout << i + 1 << " Row's maximum is " << maxarray[i] + 1 << " element with " << array[i*Cols+maxarray[i]] << " value" << endl;
		times = MPI_Wtime();
		cout << "time: " << (times - timef) << endl;
        delete[] array;
    }
    else
    {
        if (ProcRank == 0)
        {
			double timef, times, tickt;
			timef = MPI_Wtime();
			tickt = MPI_Wtick();
			int Rows, Cols;
			if (argc != 3)
				return 1;
			Rows = atoi(argv[2]);  
			Cols = atoi(argv[1]);
			if (Cols < 1 || Rows < 1)
				return 1;
			double* array = new double[Rows*Cols];   
			RandomizeArray(array, Rows*Cols);

            int RowPerProc = Rows / ProcNum;
            int MainProcRows = RowPerProc + (Rows % ProcNum);
			int* maxarray = new int[Rows];
            MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&Cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 1; i < ProcNum; i++)
				MPI_Send(array + Cols*(i-1)*RowPerProc, Cols*RowPerProc	, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MatrixFindMaxPerRow(&(array[Cols*(ProcNum-1)*RowPerProc]), Rows - (ProcNum-1)*RowPerProc, Cols, maxarray + (ProcNum - 1)*RowPerProc);
			for (int i = 1; i < ProcNum; i++)
			{
				MPI_Recv(maxarray + (i - 1)*RowPerProc, RowPerProc, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			for (int i = 0; i < Rows; i++)
			{
				cout << i + 1 << " Row's maximum is " << maxarray[i]+1 << " element with " << array[i*Cols+maxarray[i]] << " value" << endl;
			}
			times = MPI_Wtime();
			cout << "time: " << (times - timef) << endl;
            delete[] array;
			delete[] maxarray;
        }
        else
        {
            int RowPerProc;
			int Cols;
            MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&Cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
			int* maxsubarray = new int[RowPerProc];
            if (RowPerProc > 0)
            {
                double* subarray = new double[RowPerProc*Cols];
				MPI_Recv(subarray, Cols*RowPerProc, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MatrixFindMaxPerRow(subarray, RowPerProc, Cols, maxsubarray);
				MPI_Send(maxsubarray, RowPerProc, MPI_INT, 0, 0, MPI_COMM_WORLD);
                delete[] subarray;
				delete[] maxsubarray;
            }
        }

    }

    MPI_Finalize();

    return 0;
}



