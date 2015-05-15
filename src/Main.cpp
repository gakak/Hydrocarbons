#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

/*GLOBAL*/
double Psi, k, m;
double* x; int sizeX;
double* y; int sizeY;
double* z; int sizeZ;
/*GLOBAL*/

//нижний предел, верхний предел, размеры матрицы;
double** getRandMat(int min, int max, int length[2]){
	double random;
	double** p = new double* [length[0]];
	for(int i = 0; i < length[0]; i++)
		p[i] = new double [length[1]];
	for(int i = 0; i < length[0]; i++){
		random = min + rand() % max;
		for(int j = 0; j < 3; j++)
			p[i][j] = random/sqrt(3);
	}
	return p;
}

double getF(double deltaL, double g1, double c1, double c2){
	return -Psi*k*(deltaL + g1)/(c2 - c1);
}

double getS(double deltaT, double F1, double F2, double c1, double c2){
	return deltaT*(F2 - F1)/(m*(c2 - c1));
}

double** fullF(double* g_, double** L, double** F){
	//счет F
	int currentStep;
	for(int r = 0; r < sizeX; r++){
		for(int w = 0; w < sizeY; w++){
			//#pragma omp parallel for
			for(int d = 0; d < sizeZ; d++){
				currentStep = r*(sizeY*sizeZ)+w*sizeZ+d;
				if(r != 0 && r != sizeX-1)
					F[currentStep][0] = getF(L[currentStep+1][0] - L[currentStep-1][0], g_[0]*x[r], x[r-1], x[r+1]);
				else if(r == 0)
					F[currentStep][0] = getF(L[currentStep+1][0] - L[currentStep][0], g_[0]*x[r], x[r], x[r+1]);
				else if(r == sizeX-1)
					F[currentStep][0] = getF(L[currentStep][0] - L[currentStep-1][0], g_[0]*x[r], x[r-1], x[r]);

				if(w != 0 && w != sizeY-1)
					F[currentStep][1] = getF(L[currentStep+1][1] - L[currentStep-1][1], g_[1]*y[w], y[w-1], y[w+1]);
				else if(w == 0)
					F[currentStep][1] = getF(L[currentStep+1][1] - L[currentStep][1], g_[1]*y[w], y[w], y[w+1]);
				else if(w == sizeY-1)
					F[currentStep][1] = getF(L[currentStep][1] - L[currentStep-1][1], g_[1]*y[w], y[w-1], y[w]);

				if(d != 0 && d != sizeZ-1)
					F[currentStep][2] = getF(L[currentStep+1][2] - L[currentStep-1][2], g_[2]*z[d], z[d-1], z[d+1]);
				else if(d == 0)
					F[currentStep][2] = getF(L[currentStep+1][2] - L[currentStep][2], g_[2]*z[d], z[d], z[d+1]);
				else if(d == sizeZ-1)
					F[currentStep][2] = getF(L[currentStep][2] - L[currentStep-1][2], g_[2]*z[d], z[d-1], z[d]);
			}
		}
	}
	return F;
}

double** fullS(double deltaT, double** F, double** S){
	//счет S
	int currentStep;
	for(int r = 0; r < sizeX; r++){
		for(int w = 0; w < sizeY; w++){
			//#pragma omp parallel for
			for(int d = 0; d < sizeZ; d++){
				currentStep = r*(sizeY*sizeZ)+w*sizeZ+d;
				if(r != 0 && r != sizeX-1)
					S[currentStep][0] += getS(deltaT, F[currentStep-1][0], F[currentStep+1][0], x[r-1], x[r+1]);
				else if(r == 0)
					S[currentStep][0] += getS(deltaT, F[currentStep][0], F[currentStep+1][0], x[r], x[r+1]);
				else if(r == sizeX - 1)
					S[currentStep][0] += getS(deltaT, F[currentStep-1][0], F[currentStep][0], x[r-1], x[r]);

				if(w != 0 && w != sizeY-1)
					S[currentStep][1] += getS(deltaT, F[currentStep-1][1], F[currentStep+1][1], y[w-1], y[w+1]);
				else if(w == 0)
					S[currentStep][1] += getS(deltaT, F[currentStep][1], F[currentStep+1][1], y[w], y[w+1]);
				else if(w == sizeY - 1)
					S[currentStep][1] += getS(deltaT, F[currentStep-1][1], F[currentStep][1], y[w-1], y[w]);

				if(d != 0 && d != sizeZ-1)
					S[currentStep][2] += getS(deltaT, F[currentStep-1][2], F[currentStep+1][2], z[d-1], z[d+1]);
				else if(d == 0)
					S[currentStep][2] += getS(deltaT, F[currentStep][2], F[currentStep+1][2], z[d], z[d+1]);
				else if(d == sizeZ - 1)
					S[currentStep][2] += getS(deltaT, F[currentStep-1][2], F[currentStep][2], z[d-1], z[d]);
			}
		}
	}
	return S;
}

int main() {

	/*parameters in Si*/
	double t0 = 0.0, tk = 10.0, deltaT = 0.1;		// начальное, конечное время и шаг по сетке времени
	int timeSize = int((tk - t0)/deltaT + 1);		// количество точек в сетке времени
	double x0 = 1.0, xk = 10.0, deltaX = 0.1;		// начальное, конечное время и шаг по сетке координат X
	sizeX = int((xk - x0)/deltaX + 1);				// количество точек в сетке координат X
	double y0 = 1.0, yk = 5.0, deltaY = 0.1;		// начальное, конечное время и шаг по сетке координат Y
	sizeY = int((yk - y0)/deltaY + 1);				// количество точек в сетке координат Y
	double z0 = 1.0, zk = 5.5, deltaZ = 0.1;		// начальное, конечное время и шаг по сетке координат Z
	sizeZ = int((zk - z0)/deltaZ + 1);				// количество точек в сетке координат Z
	int matrixSize = sizeX*sizeY*sizeZ;				// количество точек в рассматриваемом объеме

	double** S = new double* [matrixSize];			// насыщенность(зависит от координат и времени)
	for(int i = 0; i < matrixSize; i++)
		S[i] = new double [3];

	double** F = new double* [matrixSize];			// коэффициент для вычисления насыщенности
	for(int i = 0; i < matrixSize; i++)
		F[i] = new double [3];

	m = 0.2; 										// пористость(в общем случае зависит от координат) [%/100]
	k = 100.0*pow(10, -9); 							// тензор абсолютной проницаемости(в общем случае зависит от
	//координат) [mD = 10^(-9) m^2)]
	/*petroleum*/
	double f1 = 40.55*pow(10, -9); 					// относительная фазовая проницаемость(в общем случае зависит
	//от координат) [mD = 10^(-9) m^2]
	double mu1 = 0.052; 							// динамическая вязкость [Pa*sec]
	double** p1; 									// давления(зависит от координат и времени в общем случае) [Pa]
	double ro1 = 850.0; 							// density [kg/m³]
	/*petroleum*/
	/*water*/
	double f2 = 3.75*pow(10, -9); 					// относительная фазовая проницаемость(в общем случае зависит
	//от координат) [mD = 10^(-9) m^2]
	double mu2 = 8.94*pow(10, -4); 					// динамическая вязкость [Pa*sec]
	double** p2; 									// давления(зависит от координат и времени в общем случае) [Pa]
	double ro2 = 1000.0; 							// density [kg/m³]
	/*water*/
	Psi = f1*f2/(mu1*f1 + mu2*f2); 					// коэффициент в уравнение, описывающем насыщенность одной из фаз
	double** L = new double* [matrixSize];			// разность давления в фазах
	for(int i = 0; i < matrixSize; i++)
		L[i] = new double [3];
	double g[3] = {0.0, 0.0, -9.8}; 				// вектор, описывающий объемные силы
	/*parameters*/

	/*variables*/
	x = new double [sizeX];							// grid for x
	for(int i = 0; i < sizeX; i++)
		x[i] = x0 + i*deltaX;
	y = new double [sizeY];							// grid for y
	for(int i = 0; i < sizeY; i++)
		y[i] = y0 + i*deltaY;
	z = new double [sizeZ];							// grid for z
	for(int i = 0; i < sizeZ; i++)
		z[i] = z0 + i*deltaZ;
	double* t = new double[timeSize];				// grid for time
	for(int i = 0; i < timeSize; i++)
		t[i] = t0 + i*deltaT;
	/*variables*/

	double g_[3];
	for(int i = 0; i < 3; i++){
		g_[i] = (ro2 - ro1)*g[i];
	}

	for(int i = 0; i < matrixSize; i++){
		for(int j = 0; j < 3; j++){
			F[i][j] = 1.0;
			S[i][j] = 0.0;
		}
	}


	double progtimeSt, progTime;
	int numThreads = 2, length[2] = {matrixSize, 3};
	for(int threads = 1; threads < numThreads + 1; threads++){
		omp_set_num_threads(threads);

		progtimeSt = omp_get_wtime();

		/*!!Основные расчеты!!*/
		for(int i = 0; i < timeSize-1; i++){
			/*меняем давление на каждом шаге времени*/
			p1 = getRandMat(40*98066.5, 50*98066.5, length); 					// на глубине 200m
			p2 = getRandMat(40*98066.5, 50*98066.5, length); 					// на глубине 200m
			for(int indx = 0; indx < matrixSize; indx++){
				//#pragma omp parallel for
				for(int j = 0; j < 3; j++)
					L[indx][j] = p1[indx][j] - p2[indx][j];
			}
			/*меняем давление на каждом шаге времени*/

			F = fullF(g_, L, F);
			S = fullS(t[i+1] - t[i], F, S);
		}
		/*!!Основные расчеты!!*/

		progTime = omp_get_wtime() - progtimeSt;

		/*Outputting results to the files*/
		ofstream timeFile; 														//time output file
		if(threads == 1)
			timeFile.open("time", ios_base::out | ios_base::trunc);
		else
			timeFile.open("time", ios_base::app);
		timeFile << fixed << "Время выполнения программы для " << threads << " потока(-ов): " << progTime << " sec\n";
		timeFile.close();

		ofstream file;
		file.open("saturation", ios::trunc);									//насыщенность
		if (!file.is_open())
			cout << "Can't open/create file 'saturation'!\n";
		else{
			for(int i = 0; i < matrixSize; i++){
				file << fixed << S[i][0] << "\t" << S[i][1] << "\t" << S[i][2] << "\t" << endl;
			}
		}
		file.close();

		file.open("F4saturation", ios::trunc);									//F в насыщенности
		if (!file.is_open())
			cout << "Can't open/create file 'F4saturation'!\n";
		else{
			for(int i = 0; i < matrixSize; i++){
				file << fixed << F[i][0] << "\t" << F[i][1] << "\t" << F[i][2] << "\t" << endl;
			}
		}
		file.close();
		/*Outputting results to the files*/
	}

	delete[] F;
	delete[] S;

	return 0;
}
