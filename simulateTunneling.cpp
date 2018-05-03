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
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DECLARE FUNCTIONS
// ===================================

cdouble printMatrixC(matrix k); 							// print complex matrix
double  printMatrixR(vector<double> k); 					// print real matrix
void 	saveFileC(vector<cdouble> vec, string save_name);	// save complex matrix
void 	saveFileR(vector<double> vec,  string save_name);	// save real matrix

vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);


// ===================================
// MAIN
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


