/*
 *	Created by Rachel Wooten on 7/3/2018
 *  
 *
 *
 */
#ifndef SPARSECLASS_H
#define SPARSECLASS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "MatrixClass.h"

using namespace std;

class SparseClass: public MatrixClass{
	public:

		SparseClass(){}										// Constructor
		SparseClass(const SparseClass & sp);				// Copy constructor
		SparseClass & operator= (const SparseClass & sp);  	// Copy assignment constructor
		~SparseClass();										// Destructor	

		// Symmetric matrix load and view
		void load_matrix(string infile); 	// virtual function
		void print_mat();			 		// virtual function
		void show_compression();
		double get_elem_ij(int i, int j);   // virtual function	

		// diagonalization
		void set_nev(int num_evals);
		void diagonalize();
		void dsaupd(int n, int nev, double *Evals);
		void dsaupd(int n, int nev, double *Evals, double ** Evecs);
		void av(int n, double *in, double *out);	
		void save_eigen(int how_many);

	private:
		// Store lower left triangle of matrix
		int num_ev;                	// Number of eigenvalues to be calculated; 
									// nev < dim
		
		int num_non_zero;
		double * diagonal;
		vector<double> values;		// non-zero matrix elements of lower triangle, (row-major order)
		int* first_in_row;		    // index in values of first element of row i, i = 0, 1,..., dim 
		vector<int> column_index;	// Column-index of each element in values.
		
		double zero_tol = 1e-16;	// All values x in  -zero_tol < x < zero_tol -->  0
};

#endif

