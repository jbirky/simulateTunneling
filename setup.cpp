#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

int		D 	= 601; 				// number of grid points
int 	N 	= 1; 				// number of time steps
double	h	= 1;
double	m 	= 1;
double	w	= 1;
double	X0	= -4;				// min x range
double	XD	= 4;				// max s range
double 	PI  = 3.14159265358;	
double	T0 	= 2 * PI;
int 	P 	= 32;
double	T 	= T0/16;						// period
double	DEL_T 	= T0 / 128;				// time step size
double	DEL_X	= (XD - X0)/(D - 1);	// spatial step size

// Double well potential parameters
double	ALPHA	= 1;
double 	BETA 	= 2;