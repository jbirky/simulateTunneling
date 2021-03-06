
// ===================================
// MATRIX/VECTOR OPERATION FUNCTIONS
// ===================================

cdouble printMatrixC(matrix k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


double printMatrixR(vector<double> k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


void saveFileC(vector<cdouble> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
}


void saveFileR(vector<double> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
}


vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = D;

    //Create a length D^2 array of doubles, filled with the value 0.
    vector<cdouble> matrixnew(dim*dim,0);

    //cblas zgemm takes in three matrices: A,B,C. It stores in the value C
    //the matrix alpha*AB+beta*C. In this case beta=0 and alpha=1, so we just
    //do matrix multiplication AB.

    cdouble mult(1, 0.0);
    cdouble czero(0.0, 0.0);

    //use complex matrix matrix multiply.
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &mult, &m1[0], D, &m2[0], dim, &czero, &matrixnew[0], dim);

    //return the vector object.
    return matrixnew;
}


vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = m2.size();

	vector<cdouble> prod(dim,0);

    for (int i=0; i<dim; i++) {
    	prod[i] = 0;
    	for (int j=0; j<dim; j++) {
    		prod[i] += m1[i*dim + j] * m2[j];
    	}
    }

    return prod;
}


vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = m1.size();

	vector<cdouble> prod(dim,0);

    for (int i=0; i<dim; i++) {
    	prod[i] += m1[i] * m2[i];
    }

    return prod;
}


cdouble trace(matrix K) {

    cdouble tr_sum;

    for (int i=0; i<D; i++) {
        for (int j=0; j<D; j++) {
            if (i==j) {
                tr_sum += K[i*D + j];
            }
        }
    }

    return tr_sum;
}