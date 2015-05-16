/*
 * Name			: Hydrocarbons.cpp
 * Author		: Arseniy Gerasimov, Ekaterina Mishinkina
 * Date			: 05/2015
 *
 * Theor basis	: Журавлев А. С., Вестник Тюменского государственного университета. Выпуск №7, 2011. 39-45 с.
 * 				  http://vk.cc/3N76AO
 */

#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

/*GLOBAL variables*/
double Psi, k, m;
double* x; int sizeX;
double* y; int sizeY;
double* z; int sizeZ;
/*GLOBAL variables*/

//-------------------------------------------------------------------------------------------------------------------------
double* makeGrid(double a0, double deltaA, int colsCount){
	double* g = new double [colsCount];
	for(int i = 0; i < colsCount; i++)
		g[i] = a0 + i*deltaA;
	return g;
}

//-------------------------------------------------------------------------------------------------------------------------
double** init2XArray(int rowsCount, int colsCount){
	double** a = new double* [rowsCount];
	for(int i = 0; i < rowsCount; i++)
		a[i] = new double [colsCount];
	return a;
}

//-------------------------------------------------------------------------------------------------------------------------
double** getRand2XArray(int min, int max, int rowsCount, int colsCount){
	double randNum;
	double** randArr = init2XArray(rowsCount, colsCount);
	for(int i = 0; i < rowsCount; i++){
		randNum = min + rand() % max;
		for(int j = 0; j < colsCount; j++)
			randArr[i][j] = randNum/sqrt(colsCount);
	}
	return randArr;
}

//-------------------------------------------------------------------------------------------------------------------------
double calcF(double deltaL, double g1, double c1, double c2){
	return -Psi*k*(deltaL + g1)/(c2 - c1);
}

//-------------------------------------------------------------------------------------------------------------------------
double calcS(double deltaT, double F1, double F2, double c1, double c2){
	return deltaT*(F2 - F1)/(m*(c2 - c1));
}

// calculating matrix of F ------------------------------------------------------------------------------------------------
double** getFullF(double* g_, double** L, double** F){
	int curStep;
	for(int r = 0; r < sizeX; r++){
		for(int w = 0; w < sizeY; w++){
			//#pragma omp parallel for
			for(int d = 0; d < sizeZ; d++){
				curStep = r*(sizeY*sizeZ)+w*sizeZ+d;
				if(r != 0 && r != sizeX-1)
					F[curStep][0] = calcF(L[curStep+1][0] - L[curStep-1][0], g_[0]*x[r], x[r-1], x[r+1]);
				else if(r == 0)
					F[curStep][0] = calcF(L[curStep+1][0] - L[curStep][0], g_[0]*x[r], x[r], x[r+1]);
				else if(r == sizeX-1)
					F[curStep][0] = calcF(L[curStep][0] - L[curStep-1][0], g_[0]*x[r], x[r-1], x[r]);

				if(w != 0 && w != sizeY-1)
					F[curStep][1] = calcF(L[curStep+1][1] - L[curStep-1][1], g_[1]*y[w], y[w-1], y[w+1]);
				else if(w == 0)
					F[curStep][1] = calcF(L[curStep+1][1] - L[curStep][1], g_[1]*y[w], y[w], y[w+1]);
				else if(w == sizeY-1)
					F[curStep][1] = calcF(L[curStep][1] - L[curStep-1][1], g_[1]*y[w], y[w-1], y[w]);

				if(d != 0 && d != sizeZ-1)
					F[curStep][2] = calcF(L[curStep+1][2] - L[curStep-1][2], g_[2]*z[d], z[d-1], z[d+1]);
				else if(d == 0)
					F[curStep][2] = calcF(L[curStep+1][2] - L[curStep][2], g_[2]*z[d], z[d], z[d+1]);
				else if(d == sizeZ-1)
					F[curStep][2] = calcF(L[curStep][2] - L[curStep-1][2], g_[2]*z[d], z[d-1], z[d]);
			}
		}
	}
	return F;
}

// calculating matrix of S ------------------------------------------------------------------------------------------------
double** getFullS(double deltaT, double** F, double** S){
	int curStep; //current step/текущая итерация
	for(int r = 0; r < sizeX; r++){
		for(int w = 0; w < sizeY; w++){
			//#pragma omp parallel for
			for(int d = 0; d < sizeZ; d++){
				curStep = r*(sizeY*sizeZ)+w*sizeZ+d;
				if(r != 0 && r != sizeX-1)
					S[curStep][0] += calcS(deltaT, F[curStep-1][0], F[curStep+1][0], x[r-1], x[r+1]);
				else if(r == 0)
					S[curStep][0] += calcS(deltaT, F[curStep][0], F[curStep+1][0], x[r], x[r+1]);
				else if(r == sizeX - 1)
					S[curStep][0] += calcS(deltaT, F[curStep-1][0], F[curStep][0], x[r-1], x[r]);

				if(w != 0 && w != sizeY-1)
					S[curStep][1] += calcS(deltaT, F[curStep-1][1], F[curStep+1][1], y[w-1], y[w+1]);
				else if(w == 0)
					S[curStep][1] += calcS(deltaT, F[curStep][1], F[curStep+1][1], y[w], y[w+1]);
				else if(w == sizeY - 1)
					S[curStep][1] += calcS(deltaT, F[curStep-1][1], F[curStep][1], y[w-1], y[w]);

				if(d != 0 && d != sizeZ-1)
					S[curStep][2] += calcS(deltaT, F[curStep-1][2], F[curStep+1][2], z[d-1], z[d+1]);
				else if(d == 0)
					S[curStep][2] += calcS(deltaT, F[curStep][2], F[curStep+1][2], z[d], z[d+1]);
				else if(d == sizeZ - 1)
					S[curStep][2] += calcS(deltaT, F[curStep-1][2], F[curStep][2], z[d-1], z[d]);
			}
		}
	}
	return S;
}

//-------------------------------------------------------------------------------------------------------------------------
int write3ColArray(string fileName, double** a, int rowsCount){
	ofstream outFile;
	outFile.open(fileName.c_str());
	for(int i = 0; i < rowsCount; i++){
		outFile << a[i][0] << "\t" << a[i][1] << "\t" << a[i][2] << endl;
	}
	outFile.close();
	return 1;
}

//-------------------------------------------------------------------------------------------------------------------------


int main() {

	/*parameters, given in Si*/
	double t0 = 0.0, tk = 100.0, deltaT = 0.1;	// start, end time and grid step for time
	double x0 = 1.0, xk = 3.0, deltaX = 0.1;	// start, end X and grid step for X
	double y0 = 1.0, yk = 3.0, deltaY = 0.1;	// start, end Y and grid step for Y
	double z0 = 1.0, zk = 3.0, deltaZ = 0.1;	// start, end Z and grid step for Z

	int timeSize = int((tk - t0)/deltaT + 1);	// number of points in the time grid
	sizeX = int((xk - x0)/deltaX + 1);			// number of points in the X grid
	sizeY = int((yk - y0)/deltaY + 1);			// number of points in the Y grid
	sizeZ = int((zk - z0)/deltaZ + 1);			// number of points in the Z grid
	int matrixSize = sizeX*sizeY*sizeZ;			// number of points in the grid of measured volume

	m = 0.2; 									// porosity/пористость(в общем случае зависит от координат) [%/100]
	k = 100.0*pow(10, -9); 						// absolute permeability tensor/тензор абсолютной проницаемости
	//(в общем случае зависит от координат) [mD = 10^(-9) m²)]
	double g[3] = {0.0, 0.0, -9.8}; 			// vector, describing the volume forces/вектор, описывающий объемные силы

	/*petroleum*/
	double fP = 40.55*pow(10, -9); 				// relative permeability/относительная фазовая проницаемость
	//(в общем случае зависит от координат) [mD = 10^(-9) m²]
	double muP = 0.052; 						// dynamic viscosity/динамическая вязкость [Pa*sec]
	double** pP; 								// pressure/давления(зависит от координат и времени в общем случае) [Pa]
	double roP = 850.0; 						// density/плотность [kg/m³]
	/*petroleum*/

	/*water*/
	double fW = 3.75*pow(10, -9); 				// relative permeability/относительная фазовая проницаемость
	//(в общем случае зависит от координат) [mD = 10^(-9) m²]
	double muW = 8.94*pow(10, -4); 				// dynamic viscosity/динамическая вязкость [Pa*sec]
	double** pW; 								// pressure/давления(зависит от координат и времени в общем случае) [Pa]
	double roW = 1000.0; 						// density/плотность [kg/m³]
	/*water*/

	Psi = fP*fW/(muP*fP + muW*fW); 				// коэффициент в уравнение, описывающем насыщенность одной из фаз

	double** L = init2XArray(matrixSize, 3);	// pressure difference in the phases/разность давления в фазах
	double** S = init2XArray(matrixSize, 3);	// saturation/насыщенность(зависит от координат и времени)
	double** F = init2XArray(matrixSize, 3);	// коэффициент для вычисления насыщенности
	/*parameters*/

	/*grids*/
	x = makeGrid(x0, deltaX, sizeX);			// grid for x
	y = makeGrid(y0, deltaY, sizeY);			// grid for y
	z = makeGrid(z0, deltaZ, sizeZ);			// grid for z
	double* t = makeGrid(t0, deltaT, timeSize);	// grid for time
	/*grids*/

	double g_[3];
	for(int i = 0; i < 3; i++)
		g_[i] = (roW - roP)*g[i];
	for(int i = 0; i < matrixSize; i++){
		for(int j = 0; j < 3; j++){
			F[i][j] = 1.0;
			S[i][j] = 0.0;
		}
	}


	double timeStart, runtime;
	int numThreads = 2;
	for(int threads = 1; threads < numThreads + 1; threads++){
		omp_set_num_threads(threads);

		/* Main calculations */
		timeStart = omp_get_wtime();
		for(int tIndx = 0; tIndx < timeSize-1; tIndx++){
			/*changing the pressure every time step*/
			pP = getRand2XArray(40*98066.5, 50*98066.5, matrixSize, 3);			// 200m below ground/на глубине 200m
			pW = getRand2XArray(40*98066.5, 50*98066.5, matrixSize, 3);			// 200m below ground/на глубине 200m
			for(int i = 0; i < matrixSize; i++){
				//#pragma omp parallel for
				for(int j = 0; j < 3; j++)
					L[i][j] = pP[i][j] - pW[i][j];
			}
			/*changing the pressure every time step*/

			F = getFullF(g_, L, F);
			S = getFullS(t[tIndx+1] - t[tIndx], F, S);
		}
		runtime = omp_get_wtime() - timeStart;
		/* Main calculations */

		/*Outputting results to the files*/
		ofstream timeFile;
		if(threads == 1)
			timeFile.open("time", ios_base::out | ios_base::trunc);
		else
			timeFile.open("time", ios_base::app);
		timeFile << fixed << "Время выполнения программы для " << threads << " потока(-ов): " << runtime << " sec\n";
		timeFile.close();

		if(threads == numThreads-1){
			int writeRes = write3ColArray("f4Saturation", F, matrixSize);
			if (!writeRes)
				cout << "Can't write to the file 'f4Saturation'!\n";

			writeRes = write3ColArray("saturation", S, matrixSize);
			if (!writeRes)
				cout << "Can't write to the file 'saturation'!\n";
		}
		/*Outputting results to the files*/
	}

	delete[] F;
	delete[] S;

	return 0;
}
