// ===================================
// JESSICA BIRKY (A13002163)
// ===================================

#include "setup.cpp"
#include "utils.cpp"
using namespace std;


// ===================================
// DECLARE FUNCTIONS
// ===================================

double 	x(int i);
double 	potential(double x);

// Path Integral Functions
matrix 			returnKep();
vector<cdouble> phiInit();
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

	// matrix K = returnPropagator();
	// printMatrixC(K);
	// string save_K = "prop_matrix/Kprop.dat";
	// saveFileC(K, save_K);

	// matrix phi0 = normalizeWF(phiInit());
	// vector<double> phi0_sq = returnProbability(phi0);
	// string save_phi0_prob = "wave_prob/phi_sq0.dat";
	// saveFileR(phi0_sq, save_phi0_prob);

	saveWaveFunctions();


	// matrix Kep = returnKep();
	// cout<<Kep.size()<<endl;
	// matrix K = matrixMatrixMultiply(Kep, Kep);
	// printMatrixC(K);

	// =========================================

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / pow(10,6);
	cout<<"Time elapsed: "<<elapsed_secs<<" sec"<<endl;

	return 0;
}


// ==========================================
// DEFINE POTENTIAL AND INITIAL CONDITION
// ==========================================

double x(int i)
{
	return X0 + i*DEL_X;
}


vector<cdouble> phiInit()
{
	vector<cdouble> phi0(D,0);

	for (int i=0; i<D; i++) {
		// phi0[i]  = 1 / (P_ALPHA * sqrt(2*PI)) * ( exp(-pow((x(i) - BETA),2) / (2.0*pow(P_ALPHA,2))) + exp(-pow((x(i) + BETA),2) / (2.0*pow(P_ALPHA,2))) );
		// phi0[i] += 1 / (P_ALPHA * sqrt(2*PI)) * ( exp(-pow((x(i) - BETA),2) / (2.0*pow(P_ALPHA,2))) - exp(-pow((x(i) + BETA),2) / (2.0*pow(P_ALPHA,2))) ); 
		phi0[i] = exp(-m*pow(w,2) / (2*m) * pow((x(i)-XMIN), 2));
	}

	return phi0;
}


double potential(double x)
{
	// double V = V_ALPHA * pow(x, 4) + BETA * pow(x, 2) + pow(BETA, 2) / (4*V_ALPHA);

	double V = V_ALPHA * exp(pow((pow(x,2) - pow(XMIN,2)), 2));

	return V;
}


// ===================================
// CREATE PROPAGATOR MATRIX
// ===================================

matrix returnKep()
{
	matrix Kep(D*D,0);

	// cdouble A(sqrt(PI * h * DEL_T), sqrt(PI * h * DEL_T)); 
	double A = 1; // normalize wave function later
	cdouble I(0,1);
	cdouble efactor;

	// Generate elementary Kep matrix
	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {

			efactor = I*DEL_T/h * (.5*m*pow((x(j)-x(i))/DEL_T, 2) - potential(((x(j)+x(i))/2.)));

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
		avg_pot += (conj(wf[i]) * potential(x(i)) * wf[i]).real() * DEL_X;
	}

	return avg_pot;
}
