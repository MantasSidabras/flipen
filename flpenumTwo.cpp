#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

using namespace std;

int numDP = 1000; // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;   // Esanciu objektu skaicius (preexisting facilities)
int numCL = 26;   // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 2;	 // Nauju objektu skaicius

int world_size;
int world_rank;
int name_len;

double **demandPoints; // Geografiniai duomenys

//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double *a, double *b);
double evaluateSolution(int *X, int index);
int increaseX(int *X, int index, int maxindex);
int calculateGroupCount();
void CopyArray(int *newArray, int index, int *arrayToCopy);

//=============================================================================

int main(int argc, char *argv[])
{
	double totalTime = 0;
	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the number of processes

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(processor_name, &name_len);

	for (int i = 0; i < 3; i++)
	{
		numX++;

		// printf("Hello world from processor %s, rank %d out of %d processors\n",
		// 	   processor_name, world_rank, world_size);
		double ts;
		if (world_rank == 0)
		{
			ts = getTime(); // Algoritmo vykdymo pradzios laikas
		}
		loadDemandPoints(); // Nuskaitomi duomenys
		int groupCount = calculateGroupCount();
		int leftAafterDivision = groupCount % world_size;
		int fakeGroupsCount = 0;
		if (leftAafterDivision != 0)
		{
			fakeGroupsCount = world_size - leftAafterDivision;
		}
		// if groupCount doesn't devide by thred count, add empty groups
		groupCount += fakeGroupsCount;

		int elementsPerProcessor = groupCount / world_size;
		int *allPossibleGroups;
		if (world_rank == 0)
		{
			allPossibleGroups = new int[groupCount * numX];
			for (int i = 0; i < fakeGroupsCount; i++)
			{
				//allPossibleGroups[groupCount-1-i] = new int[numX];
				for (int j = 0; j < numX; j++)
				{
					allPossibleGroups[(groupCount - 1 - i) * numX + j] = j;
				}
			}

			int *X = new int[numX];
			int *bestX = new int[numX];
			for (int i = 0; i < numX; i++)
			{
				X[i] = i;
				bestX[i] = i;
			}
			CopyArray(allPossibleGroups, 0, X);
			int groupNr = 1;

			//----- Renkamos visos grupes ------------------------------------------------

			while (true)
			{
				if (increaseX(X, numX - 1, numCL))
				{
					CopyArray(allPossibleGroups, groupNr, X);
					groupNr++;
				}
				else
					break;
			}
		}
		// allPossibleGroups << all groups
		int *subPosibleGroups = new int[elementsPerProcessor * numX];
		MPI_Scatter(allPossibleGroups, elementsPerProcessor, MPI_INT, subPosibleGroups,
					elementsPerProcessor, MPI_INT, 0, MPI_COMM_WORLD);

		double u = evaluateSolution(subPosibleGroups, 0);
		double bestU = u;
		int *bestX = new int[numX];
		for (int i = 0; i < numX; i++)
		{
			bestX[i] = subPosibleGroups[i];
		}

		for (int i = 0; i < elementsPerProcessor; i++)
		{
			u = evaluateSolution(subPosibleGroups, i);
			if (u > bestU)
			{
				bestU = u;
				for (int j = 0; j < numX; j++)
					bestX[j] = subPosibleGroups[(i * numX) + j];
			}
		}

		int *allBestX;

		if (world_rank == 0)
		{
			allBestX = new int[numX * world_size];
		}

		MPI_Gather(bestX, numX, MPI_INT, allBestX, numX, MPI_INT, 0,
				   MPI_COMM_WORLD);

		//----- Rezultatu spausdinimas --------------------------------------------
		if (world_rank == 0)
		{
			// cout << "\nThe result will blow your mind\n";
			// for (int i = 0; i < numX * world_size; i++)
			// {
			// 	cout << allBestX[i] << " ";
			// }

			int bestAnswerIndex = 0;
			double amount = evaluateSolution(allBestX, bestAnswerIndex);
			double bestAmount = amount;

			for (int i = 0; i < world_size; i++)
			{
				amount = evaluateSolution(allBestX, i);
				if (amount > bestAmount)
				{
					bestAmount = amount;
					bestAnswerIndex = i;
				}
			}

			double tf = getTime(); // Skaiciavimu pabaigos laikas

			cout << "Geriausias sprendinys: ";
			for (int i = 0; i < numX; i++)
				cout << allBestX[(bestAnswerIndex * numX) + i] << " ";
			cout << "(" << bestAmount << ")" << endl;
			cout << "Skaiciavimo trukme: " << tf - ts << endl;
			totalTime += (tf - ts);
		}
	}
	if (world_rank == 0)
	{
		cout << "total calc time: " << totalTime << endl;
	}
	MPI_Finalize();
}

//=============================================================================

void loadDemandPoints()
{
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double *[numDP];
	for (int i = 0; i < numDP; i++)
	{
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double *a, double *b)
{
	double dlon = fabs(a[0] - b[0]);
	double dlat = fabs(a[1] - b[1]);
	double aa = pow((sin((double)dlon / (double)2 * 0.01745)), 2) + cos(a[0] * 0.01745) * cos(b[0] * 0.01745) * pow((sin((double)dlat / (double)2 * 0.01745)), 2);
	double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
	double d = 6371 * c;
	return d;
}

//=============================================================================

double getTime()
{
	struct timeval laikas;
	gettimeofday(&laikas, NULL);
	double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
	return rez;
}

//=============================================================================

double evaluateSolution(int *X, int index)
{
	double U = 0;
	int bestPF;
	int bestX;
	double d;

	for (int i = 0; i < numDP; i++)
	{
		bestPF = 1e5;
		for (int j = 0; j < numPF; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF)
				bestPF = d;
		}
		bestX = 1e5;
		for (int j = 0; j < numX; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[X[(index * numX) + j]]);

			if (d < bestX)
				bestX = d;
		}
		if (bestX < bestPF)
			U += demandPoints[i][2];
		else if (bestX == bestPF)
			U += 0.3 * demandPoints[i][2];
	}
	return U;
}

//=============================================================================
// Jei pavyko papliusinti, tai 1, o jei ne, tai 0
int increaseX(int *X, int index, int maxindex)
{
	// Ar x mazesnis uz candidate location? x[4] < 26, x[3] < 25, x[2] < 24 ...
	// ar A < B. B visada nuo 26 iki 22
	if (X[index] + 1 < maxindex - (numX - index - 1))
	{
		X[index]++;
	}
	else
	{
		if ((index == 0) && (X[index] + 1 == maxindex - (numX - index - 1)))
		{
			return 0;
		}
		else
		{
			if (increaseX(X, index - 1, maxindex))
				X[index] = X[index - 1] + 1;
			else
				return 0;
		}
	}
	return 1;
}

int calculateGroupCount()
{
	int count = numCL;
	for (int i = numCL - 1; i > numCL - numX; i--)
	{
		count *= i;
	}

	for (int i = 1; i <= numX; i++)
	{
		count /= i;
	}
	return count;
}

void CopyArray(int *newArray, int index, int *arrayToCopy)
{
	for (int i = 0; i < numX; i++)
	{
		newArray[(index * numX) + i] = arrayToCopy[i];
	}
}
