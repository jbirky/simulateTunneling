#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

int		D 	= 601; 				// number of grid points
int 	N 	= 8; 				// number of time steps
double	h	= 1;
double	m 	= 1;
double	w	= 1;
double	XS	= 0.75;				// center of initial condition gaussian
double	X0	= -4;				// min x range
double	XD	= 4;				// max s range
double 	PI  = 3.14159265358;	
double	T0 	= 2 * PI;
int 	P 	= 16;
double	T 	= T0/P;						// period
double	DEL_T 	= T0 / 128;				// time step size
double	DEL_X	= (XD - X0)/(D - 1);	// spatial step size
double	ALPHA	= 2;

double sin(double x);
double cos(double x);


// ===================================
// DECLARE FUNCTIONS
// ===================================

vector<cdouble> phiInit();
matrix 			returnKep();
matrix 			returnPropagator();							// return propagator matrix at time step N
void			saveWaveFunctions();						// save wave functions, probability amplitudes, and expected values at each time step
vector<double> 	returnProbability(vector<cdouble> wf);		// return wave function squared
vector<cdouble> normalizeWF(vector<cdouble> wf);			// normalize wave function such that sum(wf* wf) = 1

double 	avgPosition(vector<cdouble> wf); 					// return <x> for each time step up to n
double 	avgKinetic(vector<cdouble> wf);	 					// return <K> for each time step up to n
double 	avgPotential(vector<cdouble> wf); 					// return <V> for each time step up to n


// ===================================
// CREATE PROPAGATOR MATRIX
// ===================================


matrix returnKep()
{
	matrix Kep(D*D,0);

	// cdouble A(sqrt(PI * h * DEL_T), sqrt(PI * h * DEL_T)); 
	double A = 1; // normalize wave function later
	double V = .5*m*pow(w,2);
	cdouble I(0,1);
	cdouble efactor;

	// Generate elementary Kep matrix
	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {

			efactor = I*DEL_T/h * (.5*m*pow((x(j)-x(i))/DEL_T, 2) - V * pow(((x(j)+x(i))/2.), 2));
			// cout<<exp(efactor)<<endl;

			Kep[i*D + j] = exp(efactor) / A;
		}
	}

	return Kep;
}


matrix returnPropagator()
{
	// Return the Nth step propagator matrix

	matrix Kep = returnKep();
	matrix K = Kep;

	// Multiply to the nth power to propagate by n time steps
	for (int n=1; n<N; n++) {
		K = matrixMatrixMultiply(K, Kep);
	}

	// Multiply by step size factors
	double factor = pow(DEL_X, N-1);// * DEL_T;

	for (int i=0; i<K.size(); i++) {
		K[i] *= factor;
	}

	return K;
}


// ===================================
// COMPUTE WAVE FUNCTION VECTOR
// ===================================

vector<cdouble> normalizeWF(vector<cdouble> wf)
{

	vector<cdouble> norm_wf(D,0);
	double sum = 0;

	// Square wave function: multiply by complex conjugate
	for (int i=0; i<D; i++) {
		sum += (wf[i] * conj(wf[i])).real() * DEL_X;
	}

	// Divide by normalization factor such that sum(phi*phi) = 1
	for (int i=0; i<D; i++) {
		norm_wf[i] = wf[i] / sqrt(sum);
	}

	return norm_wf;
}


vector<double> returnProbability(vector<cdouble> wf) 
{
	vector<double> wf_sq(D,0);

	for (int i=0; i<D; i++) {
		wf_sq[i] = (conj(wf[i]) * wf[i]).real();
	}

	return wf_sq;
}


void saveWaveFunctions()
{
	vector<cdouble> wf_n;
	vector<cdouble> wf_n_1;
	vector<cdouble> wf_norm;
	vector<double> wf_sq;
	string sname0;
	string sname1;
	vector<double> avg_pos(P+1,0);
	vector<double> avg_eng(P+1,0);
	vector<double> avg_kin(P+1,0);
	vector<double> avg_pot(P+1,0);

	matrix K = returnPropagator();

	for (int i=0; i<P+1; i++) {

		if (i == 0) {
			wf_n = phiInit();
			wf_n_1 = wf_n;
		} else {
			wf_n = matrixVectorMultiply(K, wf_n_1);
			wf_n_1 = wf_n;
		}

		// Save wave function vectors
		wf_norm = normalizeWF(wf_n);
		sname0 = "wave_func/phi" + to_string(i) + ".dat";
		saveFileC(wf_norm, sname0);

		// Save probability amplitudes
		wf_sq = returnProbability(wf_norm);
		sname1 = "wave_prob/phi_sq" + to_string(i) + ".dat";
		saveFileR(wf_sq, sname1);

		// Compute average values
		avg_pos[i] = avgPosition(wf_norm);
		avg_pot[i] = avgPotential(wf_norm);
		avg_kin[i] = avgKinetic(wf_norm);
		avg_eng[i] = avg_pot[i] + avg_kin[i];
	}

	// Save average values
	string sname2 = "expected/avg_pos.dat";
	string sname3 = "expected/avg_pot.dat";
	string sname4 = "expected/avg_kin.dat";
	string sname5 = "expected/avg_eng.dat";
	saveFileR(avg_pos, sname2);
	saveFileR(avg_pot, sname3);
	saveFileR(avg_kin, sname4);
	saveFileR(avg_eng, sname5);
}


// ===================================
// EXPECTED VALUES
// ===================================

double avgPosition(vector<cdouble> wf)
{
	double avg_pos = 0;
	
	for (int i=0; i<D; i++) {
		avg_pos += (conj(wf[i]) * x(i) * wf[i]).real() * DEL_X;
	}

	return avg_pos;
}


double avgKinetic(vector<cdouble> wf)
{
	double avg_kin = 0;

	for (int i=0; i<D; i++) {
		if (i == 0) {
			avg_kin += (-h*h / (2*m)) * (conj(wf[i]) * (wf[i+1] - 2.*wf[i])).real() / DEL_X;
		} else if (i == D) {
			avg_kin += (-h*h / (2*m)) * (conj(wf[i]) * (- 2.*wf[i] + wf[i-1])).real() / DEL_X;
		} else {
			avg_kin += (-h*h / (2*m)) * (conj(wf[i]) * (wf[i+1] - 2.*wf[i] + wf[i-1])).real() / DEL_X;
		}
	}

	return avg_kin;
}


double avgPotential(vector<cdouble> wf)
{
	double avg_pot = 0;

	for (int i=0; i<D; i++) {
		avg_pot += (conj(wf[i]) * .5*m*w*pow(x(i),2) * wf[i]).real() * DEL_X;
	}

	return avg_pot;
}
