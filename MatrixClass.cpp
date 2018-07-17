/*
 *	Implementation of MatrixClass
 *
 */

#include "MatrixClass.h"

void MatrixClass::find_dim(){
	int count = 0;
	string line;
	ifstream file(filename);
	while( getline(file,line))
		count++;
	dim = count;
	file.close();

	return;
}

void MatrixClass::allocate_evals(){
	eigenvalues = new double [dim];
	eigenvectors = new double * [dim];
	for( int i = 0; i < dim; i++)
		eigenvectors[i] = new double [dim];
	return;
}


