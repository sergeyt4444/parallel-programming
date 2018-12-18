#include <iostream>
#include <iomanip> 
#include "time.h"
#include <cstdlib>
#include <cmath>
#include "mpi.h"
#define step 1000
using namespace std;

void RandomizeArray(double* rarray, int size, int min = -50, int max = 50)
{
	for (int i = 0; i < size; i++)
	{
		*(rarray + i) = rand() % (max - min + 1) + min;
	}
}


void PrintPoints(double* X, double* Y, int Size)
{
	for (int i = 0; i < Size; i++)
	{
		cout << setw(3) << X[i] << ", " << Y[i] << "; ";
	}
}


void ReinitEqPoints(double* X, double* Y, int Size, int min = -50, int max = 50)
{
	for (int i = 0; i < Size - 1; i++)
		for (int j = i + 1; j < Size; j++)
		{
			if (X[i] == X[j] && Y[i] == Y[j])
			{
				int flag = 1;
				while (flag)
				{
					flag = 0;
					X[i] = rand() % (max - min + 1) + min;
					Y[i] = rand() % (max - min + 1) + min;
					for (int k = 0; k < Size; k++)
					{
						if (X[i] == X[k] && Y[i] == Y[k] && i != k)
						{
							flag = 1;
							break;
						}
					}
				}
				break;
			}
		}
}

int FindBLPoint(double* X, double* Y, int Size)
{
	if (Size < 1)
	{
		return -1;
	}
	int pointInd = 0;
	for (int i = 1; i < Size; i++)
	{
		if (Y[i] < Y[pointInd])
			pointInd = i;
		else
		{
			if (Y[i] == Y[pointInd] && X[i] < X[pointInd])
				pointInd = i;
		}
	}
	return pointInd;
}

int FindTRPoint(double* X, double* Y, int Size)
{
	if (Size < 1)
	{
		return -1;
	}
	int pointInd = 0;
	for (int i = 1; i < Size; i++)
	{
		if (Y[i] > Y[pointInd])
			pointInd = i;
		else
		{
			if (Y[i] == Y[pointInd] && X[i] > X[pointInd])
				pointInd = i;
		}
	}
	return pointInd;
}

double FindDist(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

double GetCos(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return ((x2 - x1)*(x3 - x1) + (y2 - y1)*(y3 - y1)) / (FindDist(x1, y1, x2, y2)*FindDist(x1, y1, x3, y3));
}


int FindPWithMinAngle(double* X, double* Y, int Size, double x1, double y1, double x2, double y2)
{
	double maxCos = -1.5;
	int NextPointInd = -1;
	double tmp;
	for (int i = 0; i < Size; i++)
	{
		if (!((X[i] == x1 && Y[i] == y1) || (X[i] == x2 && Y[i] == y2)))
		{
			if (((x2 - x1)*(Y[i] - y2) - (X[i] - x2)*(y2 - y1)) >= 0)
			{
				tmp = (-1)*GetCos(x2, y2, x1, y1, X[i], Y[i]);
				if (tmp > maxCos)
				{
					maxCos = tmp;
					NextPointInd = i;
				}
			}
		}
	}
	return NextPointInd;
}

int* FindEnvLinear(double* X_coord, double* Y_coord, int Size, int& PNum)
{
	int FirstPoint = FindBLPoint(X_coord, Y_coord, Size);
	int* Envelope = new int[Size + 1];
	Envelope[0] = FirstPoint;
	Envelope[1] = FindPWithMinAngle(X_coord, Y_coord, Size, X_coord[FirstPoint] - 1, Y_coord[FirstPoint], X_coord[FirstPoint], Y_coord[FirstPoint]);
	PNum = 1;
	while (Envelope[PNum] != FirstPoint)
	{
		PNum++;
		Envelope[PNum] = FindPWithMinAngle(X_coord, Y_coord, Size, X_coord[Envelope[PNum - 2]], Y_coord[Envelope[PNum - 2]],
			X_coord[Envelope[PNum - 1]], Y_coord[Envelope[PNum - 1]]);
	}
	return Envelope;
}

void ElimPointsOnLines(double* X, double* Y, int* Envelope, int& Size)
{
	int i = 1;
	while (i != Size)
	{
		if ((double)((double)(X[Envelope[i]] - X[Envelope[i-1]])*(double)(Y[Envelope[i+1]] - Y[Envelope[i]])) == (double)((double)(X[Envelope[i+1]] - X[Envelope[i]])*(double)(Y[Envelope[i]] - Y[Envelope[i-1]])))
		{
			int j = 0;
			while (i + j < Size)
			{
				Envelope[i+j] = Envelope[i+j + 1];
				j++;
			}
			Size--;
			i--;
		}
		i++;
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
		double timef, times, timepart, tickt;
		timef = MPI_Wtime();
		tickt = MPI_Wtick();
		int Size;
		if (argc != 2 && argc!=4)
			return 1;
		Size = atoi(argv[1]);
		if (Size < 1)
			return 1;
		double* X_coord = new double[Size];
		double* Y_coord = new double[Size];
		if (argc == 4)
		{
		int minRand, maxRand;
			minRand = atoi(argv[2]);
			maxRand = atoi(argv[3]);
			if (maxRand < minRand)
				return 1;
			RandomizeArray(X_coord, Size, minRand, maxRand);
			RandomizeArray(Y_coord, Size, minRand, maxRand);
			if (Size > (maxRand - minRand + 1)*(maxRand - minRand + 1))
				return 1;
			ReinitEqPoints(X_coord, Y_coord, Size, minRand, maxRand);
		}
		else 
		{
			RandomizeArray(X_coord, Size);
			RandomizeArray(Y_coord, Size);
			if (Size > 101 * 101)
				return 1;
			ReinitEqPoints(X_coord, Y_coord, Size);
		}
		PrintPoints(X_coord, Y_coord, Size);
		cout << endl;
		timepart = MPI_Wtime();
		if (Size == 1)
		{
			cout << "Result chain of points is a single point:" << endl;
			cout << setw(3) << X_coord[0] << ", " << Y_coord[0] << "; ";
		}
		else
		{
			if (Size == 2)
			{
				cout << "Result chain of points is" << endl;
				cout << setw(3) << X_coord[0] << ", " << Y_coord[0] << "; ";
				cout << setw(3) << X_coord[1] << ", " << Y_coord[1] << "; ";
			}
			else
			{
				int* Envelope = new int[Size+1];
				int PNum = 0;
				Envelope = FindEnvLinear(X_coord, Y_coord, Size, PNum);
				ElimPointsOnLines(X_coord, Y_coord,Envelope, PNum);
				if (PNum < 3)
				{
					cout << "Result chain of points is a line" << endl;
					int point;
					point = FindBLPoint(X_coord, Y_coord, Size);
					cout << setw(3) << X_coord[point] << ", " << Y_coord[point] << "; ";
					point = FindTRPoint(X_coord, Y_coord, Size);
					cout << setw(3) << X_coord[point] << ", " << Y_coord[point] << "; ";
				}
				else
				{
					cout << "Result chain of points is" << endl;
					for (int i = 0; i < PNum; i++)
					{
						cout << setw(3) << X_coord[Envelope[i]] << ", " << Y_coord[Envelope[i]] << "; ";
					}
				}
				delete[] Envelope;
			}
		}
		times = MPI_Wtime();
		cout << "time: " << (times - timef) << endl;
		cout << "time without initialisation and preparations: " << (times-timepart) << endl;
		delete[] X_coord;
		delete[] Y_coord;
	}
	else
	{
		if (ProcRank == 0)
		{

			double timef, times, timepart, timepartf, tickt;
			timef = MPI_Wtime();
			tickt = MPI_Wtick();
			int Size;
			if (argc != 2 && argc != 4)
				return 1;
			Size = atoi(argv[1]);
			if (Size < 1)
				return 1;
			double* X_coord = new double[Size];
			double* Y_coord = new double[Size];
			if (argc == 4)
			{
				int minRand, maxRand;
				minRand = atoi(argv[2]);
				maxRand = atoi(argv[3]);
				if (maxRand < minRand)
					return 1;
				RandomizeArray(X_coord, Size, minRand, maxRand);
				RandomizeArray(Y_coord, Size, minRand, maxRand);
				if (Size >(maxRand - minRand + 1)*(maxRand - minRand + 1))
					return 1;
				ReinitEqPoints(X_coord, Y_coord, Size, minRand, maxRand);
			}
			else
			{
				RandomizeArray(X_coord, Size);
				RandomizeArray(Y_coord, Size);
				if (Size > 101 * 101)
					return 1;
				ReinitEqPoints(X_coord, Y_coord, Size);
			}
			PrintPoints(X_coord, Y_coord, Size);
			timepart = MPI_Wtime();
			cout << endl;
			MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (Size == 1)
			{
				cout << "Result chain of points is a single point:" << endl;
				cout << setw(3) << X_coord[0] << ", " << Y_coord[0] << "; ";
				
			}
			else
			{
				if (Size == 2)
				{
					cout << "Result chain of points is" << endl;
					cout << setw(3) << X_coord[0] << ", " << Y_coord[0] << "; ";
					cout << setw(3) << X_coord[1] << ", " << Y_coord[1] << "; ";
				}
				else
				{
					int PointPerProc = Size / ProcNum;
					int basesize = PointPerProc + Size%ProcNum;
					MPI_Bcast(&PointPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&basesize, 1, MPI_INT, 0, MPI_COMM_WORLD);
					double* X_loc = new double[basesize];
					double* Y_loc = new double[basesize];
					double* X_buf = new double[ProcNum];
					double* Y_buf = new double[ProcNum];
					int EnvPoint;
					int* Envelope;
					int dynsize = step;
					Envelope = (int*)malloc(sizeof(int) * dynsize);
					int* buf = new int[ProcNum];
					int* sizes = new int[ProcNum];
					int* place = new int[ProcNum];
					sizes[0] = basesize;
					place[0] = 0;
					for (int i = 1; i < ProcNum; i++)
					{
						sizes[i] = PointPerProc;
						place[i] = place[i - 1] + sizes[i - 1];
					}
					MPI_Scatterv(X_coord, sizes, place, MPI_DOUBLE, X_loc, basesize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(Y_coord, sizes, place, MPI_DOUBLE, Y_loc, basesize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					int FirstPoint, LocalFP;
					LocalFP = FindBLPoint(X_loc, Y_loc, basesize);
					MPI_Gather(&LocalFP, 1, MPI_INT, buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
					for (int i = 0; i < ProcNum; i++)
					{
						if (buf[i] < 0)
						{
							buf[i] = buf[0];
						}
						X_buf[i] = X_coord[buf[i]];
						Y_buf[i] = Y_coord[buf[i]];
					}
					FirstPoint = buf[FindBLPoint(X_buf, Y_buf, ProcNum)];
					int PNum = 1;
					Envelope[0] = FirstPoint;
					MPI_Bcast(X_coord + FirstPoint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(Y_coord + FirstPoint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					EnvPoint = FindPWithMinAngle(X_loc, Y_loc, basesize, X_coord[FirstPoint]-1, Y_coord[FirstPoint], X_coord[FirstPoint] , Y_coord[FirstPoint]);
					MPI_Gather(&EnvPoint, 1, MPI_INT, buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
					for (int i = 0; i < ProcNum; i++)
					{
						if (buf[i] < 0)
						{
							buf[i] = buf[0];
						}
						X_buf[i] = X_coord[buf[i]];
						Y_buf[i] = Y_coord[buf[i]];
					}
					Envelope[1] = buf[FindPWithMinAngle(X_buf, Y_buf, ProcNum, X_coord[FirstPoint]-1, Y_coord[FirstPoint], X_coord[FirstPoint], Y_coord[FirstPoint])];
					while (Envelope[PNum] != FirstPoint)
					{
						PNum++;
						MPI_Bcast(X_coord + Envelope[PNum-1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
						MPI_Bcast(Y_coord + Envelope[PNum-1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
						EnvPoint = FindPWithMinAngle(X_loc, Y_loc, basesize, X_coord[Envelope[PNum - 2]], Y_coord[Envelope[PNum - 2]],
							X_coord[Envelope[PNum - 1]], Y_coord[Envelope[PNum - 1]]);
						MPI_Gather(&EnvPoint, 1, MPI_INT, buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
						for (int i = 0; i < ProcNum; i++)
						{
							if (buf[i] < 0)
							{
								buf[i] = buf[0];
							}
							X_buf[i] = X_coord[buf[i]];
							Y_buf[i] = Y_coord[buf[i]];
						}
						if (PNum == dynsize)
						{
							dynsize += step;
							Envelope = (int*)realloc(Envelope, sizeof(int) * (dynsize));
						}
						Envelope[PNum] = buf[FindPWithMinAngle(X_buf, Y_buf, ProcNum, X_coord[Envelope[PNum - 2]], Y_coord[Envelope[PNum - 2]],
							X_coord[Envelope[PNum - 1]], Y_coord[Envelope[PNum - 1]])];
					}
					MPI_Bcast(X_coord + Envelope[PNum - 1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(Y_coord + Envelope[PNum - 1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					ElimPointsOnLines(X_coord, Y_coord, Envelope, PNum);
					timepartf = MPI_Wtime();
					int* EnvelopeForCheck = new int[Size + 1];
					int PNumForCheck = 0;
					EnvelopeForCheck = FindEnvLinear(X_coord, Y_coord, Size, PNumForCheck);
					ElimPointsOnLines(X_coord, Y_coord, EnvelopeForCheck, PNumForCheck);
					int Correct = 1;
					if (PNum != PNumForCheck)
						Correct = 0;
					else 
					{
						for (int i = 0; i < PNum; i++)
							if (Envelope[i] != EnvelopeForCheck[i])
							{
								Correct = 0;
								break;
							}
					}
					if (Correct)
						cout << "Results checked successfully" << endl;
					else
						cout << "Results dont match" << endl;
					if (PNum < 3)
					{
						cout << "Result chain of points is a line" << endl;
						int point;
						point = FindBLPoint(X_coord, Y_coord, Size);
						cout << setw(3) << X_coord[point] << ", " << Y_coord[point] << "; ";
						point = FindTRPoint(X_coord, Y_coord, Size);
						cout << setw(3) << X_coord[point] << ", " << Y_coord[point] << "; ";
					}
					else
					{
						cout << "Result chain of points is" << endl;
						for (int i = 0; i < PNum; i++)
						{
							cout << setw(3) << X_coord[Envelope[i]] << ", " << Y_coord[Envelope[i]] << "; ";
						}
					}
					times = MPI_Wtime();
					cout << "time: " << (times - timef) << endl;
					cout << "time without initialisation and preparations: " << (timepartf-timepart) << endl;
					free(Envelope);
					delete[] EnvelopeForCheck;
					delete[] X_coord;
					delete[] X_loc;
					delete[] X_buf;
					delete[] Y_coord;
					delete[] Y_loc;
					delete[] Y_buf;
					delete[] sizes;
					delete[] place;
				}
			}
		}
		else
		{
			int Size;
			MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (Size > 2)
			{
				int PointPerProc;
				MPI_Bcast(&PointPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
				int basesize;
				MPI_Bcast(&basesize, 1, MPI_INT, 0, MPI_COMM_WORLD);
				double* X_coord = new double[PointPerProc];
				double* Y_coord = new double[PointPerProc];
				double x1, x2, y1, y2;
				int EnvPoint;
				MPI_Scatterv(0, 0, 0, MPI_DOUBLE, X_coord, PointPerProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(0, 0, 0, MPI_DOUBLE, Y_coord, PointPerProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				int LocalFP = FindBLPoint(X_coord, Y_coord, PointPerProc);
				if (LocalFP >= 0)
				{
					LocalFP += basesize + (ProcRank - 1)* PointPerProc;
				}
				else
				{
					LocalFP = -1;
				}
				MPI_Gather(&LocalFP, 1, MPI_INT, 0, 0, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&x1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&y1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				EnvPoint = FindPWithMinAngle(X_coord, Y_coord, PointPerProc, x1-1, y1, x1 , y1);
				if (EnvPoint >= 0)
				{
					EnvPoint += basesize + (ProcRank - 1)* PointPerProc;
				}
				else
				{
					EnvPoint = -1;
				}
				MPI_Gather(&EnvPoint, 1, MPI_INT, 0, 0, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&x2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&y2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				while (x1 != x2 || y1 != y2)
				{
					EnvPoint = FindPWithMinAngle(X_coord, Y_coord, PointPerProc, x1, y1, x2, y2);
					if (EnvPoint >= 0)
					{
						EnvPoint += basesize + (ProcRank - 1)* PointPerProc;
					}
					else
					{
						EnvPoint = -1;
					}
					MPI_Gather(&EnvPoint, 1, MPI_INT, 0, 0, MPI_INT, 0, MPI_COMM_WORLD);
					x1 = x2;
					y1 = y2;
					MPI_Bcast(&x2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(&y2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				}
				delete[] X_coord;
				delete[] Y_coord;
			}
		}
	}
	MPI_Finalize();
	return 0;
}