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

void PrintVector(double* arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << arr[i] << " ";
	}
}

void PrintMatrix(double* matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		cout << endl;
		PrintVector(matrix + i*cols, cols);
	}
}

double MultRowByVec(double* mrow, double* vec, int size)
{
	double res = 0;
	for (int i = 0; i < size; i++)
	{
		res += mrow[i] * vec[i];
	}
	return res;
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
		double* vector = new double[Cols];
		RandomizeArray(vector, Cols);
		PrintVector(vector, Cols);
		cout << endl;
		PrintMatrix(array, Rows, Cols);
		double* result = new double[Rows];
		for (int i = 0; i < Rows; i++)
		{
			result[i] = MultRowByVec(array + i*Cols, vector, Cols);
		}
		cout << endl << "Result vector is" << endl;
		for (int i = 0; i < Rows; i++)
			cout << result[i] << " ";
		times = MPI_Wtime();
		cout << "time: " << (times - timef) << endl;
		delete[] vector;
		delete[] result;
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
			double* vector = new double[Cols];
			RandomizeArray(vector, Cols);
			PrintVector(vector, Cols);
			cout << endl;
			PrintMatrix(array, Rows, Cols);
			int RowPerProc = Rows / ProcNum;
			double* result = new double[Rows];
			MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&Cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(vector, Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			int basesize = Cols* (RowPerProc + Rows%ProcNum);
			double* buf = new double[basesize];
			int* sizes = new int[ProcNum];
			int* place = new int[ProcNum];
			double* res = new double[RowPerProc + Rows%ProcNum];
			sizes[0] = basesize;
			place[0] = 0;
			for (int i = 1; i < ProcNum; i++)
			{
				sizes[i] = RowPerProc*Cols;
				place[i] = place[i - 1] + sizes[i - 1];
			}
			MPI_Scatterv(array, sizes, place, MPI_DOUBLE, buf, basesize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			for (int i = 0; i < RowPerProc + Rows%ProcNum; i++)
			{
				res[i] = MultRowByVec(buf + Cols*i, vector, Cols);
			}
			sizes[0] = RowPerProc + Rows%ProcNum;
			place[0] = 0;
			for (int i = 1; i < ProcNum; i++)
			{
				sizes[i] = RowPerProc;
				place[i] = place[i - 1] + sizes[i - 1];
			}
			MPI_Gatherv(res, sizes[0], MPI_DOUBLE, result, sizes, place, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			cout << endl << "Result vector is" <<endl;
			for (int i = 0; i < Rows; i++)
			{
				cout  << result[i] << "  ";
			}
			times = MPI_Wtime();
			cout << endl;
			cout << "time: " << (times - timef) << endl;
			delete[] array;
			delete[] vector;
			delete[] result;
			delete[] buf;
			delete[] sizes;
			delete[] place;
			delete[] res;
		}
		else
		{
			int RowPerProc;
			int Cols;
			MPI_Bcast(&RowPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&Cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
			double* vector = new double[Cols];
			MPI_Bcast(vector, Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (RowPerProc > 0)
			{
				double* res = new double[RowPerProc];
				double* buf = new double[RowPerProc*Cols];
				MPI_Scatterv(0, 0, 0, MPI_DOUBLE, buf, Cols*RowPerProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				for (int i = 0; i < RowPerProc; i++)
				{
					res[i] = MultRowByVec(buf + Cols*i, vector, Cols);
				}
				MPI_Gatherv(res, RowPerProc, MPI_DOUBLE, 0, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				delete[] res;
				delete[] buf;
			}
			delete[] vector;
		}

	}

	MPI_Finalize();

	return 0;
}