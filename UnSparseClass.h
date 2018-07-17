/*
 *	Created by Rachel Wooten on 7/3/2018
 *
 */

#ifndef UNSPARSECLASS_H
#define UNSPARSECLASS_H	

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
//#include <Accelerate/Accelerate.h>

#include "MatrixClass.h"

using namespace std;

class UnSparseClass: public MatrixClass{

	public:
		UnSparseClass()=default;								// Constructor
		UnSparseClass(const UnSparseClass &un);					// Copy constructor
		UnSparseClass & operator=(const UnSparseClass & un);	// Copy assignment constructor
		~UnSparseClass();										// Destructor
		


		// Symmetric matrix load and view
		void load_matrix(string infile);
		void print_mat();
		double get_elem_ij(int i, int j);


		// Diagonalization (LAPACK interface)
		void diagonalize( );
		void save_eigen(int how_many);

	protected:
		
	private:
		double * array;   // Symmetric matrix stored in Upper packed storage

};

#endif
