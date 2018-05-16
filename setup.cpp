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

// Physical constants
double	h	= 1;
double	m 	= 1;
double	w	= 1;
double 	PI  = 3.14159265358;	

// Grid parameters
int		D 	= 601; 						// number of grid points
double	X0	= -6;						// min x range
double	XD	= 6;						// max s range
double	DEL_X	= (XD - X0)/(D - 1);	// spatial step size

// Time evolution parameters
int 	N 	= 8; 						// number of time steps per K matrix
double	T0 	= 2 * PI;
int 	P 	= 200;
double	T 	= T0/16;					// period
double	DEL_T = T0 / 128;				// time step size

// Intitial phi parameters
double P_ALPHA = 1.25;
double XMIN	   = 2.5;

// Double well potential parameters
double	V_ALPHA	= .02;
double 	BETA 	= 2;