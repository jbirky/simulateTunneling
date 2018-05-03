// ===================================
// JESSICA BIRKY (A13002163)
// ===================================

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>
#include <cblas.h>
#include "utils.cpp"
#include "pathIntegral.cpp"
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DECLARE FUNCTIONS
// ===================================

// From pathIntegral.cpp
matrix 			returnKep();
matrix 			returnPropagator();							// return propagator matrix at time step N
void			saveWaveFunctions();						// save wave functions, probability amplitudes, and expected values at each time step
vector<double> 	returnProbability(vector<cdouble> wf);		// return wave function squared
vector<cdouble> normalizeWF(vector<cdouble> wf);			// normalize wave function such that sum(wf* wf) = 1

double 	avgPosition(vector<cdouble> wf); 					// return <x> for each time step up to n
double 	avgKinetic(vector<cdouble> wf);	 					// return <K> for each time step up to n
double 	avgPotential(vector<cdouble> wf); 					// return <V> for each time step up to n

// From utils.cpp
cdouble printMatrixC(matrix k); 							// print complex matrix
double  printMatrixR(vector<double> k); 					// print real matrix
void 	saveFileC(vector<cdouble> vec, string save_name);	// save complex matrix
void 	saveFileR(vector<double> vec,  string save_name);	// save real matrix

vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);


// ===================================
// RUN SIMULATION
// ===================================

int main() 
{
	clock_t begin = clock();

	// ========================================= 



	// =========================================

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / pow(10,6);
	cout<<"Time elapsed: "<<elapsed_secs<<" sec"<<endl;

	return 0;
}


// ==========================================
// DEFINE POTENTIAL AND INITIAL CONDITION
// ==========================================

double x(int d)
{
	return X0 + d*DEL_X;
}


vector<cdouble> phiInit()
{
	vector<cdouble> phi0(D,0);

	for (int i=0; i<D; i++) {
		phi0[i] = pow(ALPHA/PI, .25) * exp(-ALPHA * pow((x(i) - XS), 2) /2);
	}

	return phi0;
}


vector<cdouble> potential()
{
	vector<cdouble> V(D,0);

	return V;
}


// ===================================
// MONTE CARLO 
// ===================================