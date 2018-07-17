/*
 *  Created by Rachel Wooten, 7/3/2018
 *
 *	This code will pack a symmetric sparse matrix with a non-zero diagonal, and
 *	diagonalize it in c++ using the FORTRAN version of ARPACK.
 *
 *	-- results will be tested by comparing to the results from the 
 *	   gold standard numerical linear algebra package, LAPACK
 *  -- 
 *
 *  Classes:
 *  	MatrixClass--  General (abstract) base class to store general matrix
 *  				  properties
 *  	UnSparseClass:MatrixClass--  Specialized subclass for ordinary, 
 *  	                             uncompressed matrices
 *  	SparseClass:MatrixClass-- Specialized subclass for compressee matrices
 *
 *	IMPORTANT NOTE:
 *		In Fortran, matrices are stored in column-major order, so in standard
 *		matrix notation, element from matrix A from row i, column j is denoted
 *		A(j, i)
 *		and the "fast" index is the first index
 *
 *		This is in contrast to C++, where matrices are stored in row-major order, 
 *		which is more typical for mathematics,
 *		A[i][j]
 *		where the "fast" index is the last index.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
//#include <Accelerate/Accelerate.h>

#include "MatrixClass.h"
#include "UnSparseClass.h"
#include "SparseClass.h"

using namespace std;

int main(){

	string mat_file = "matrix.dat";
//	cout << "Dimension of matrix_in is :" << find_dim(mat_file) << endl;	

	UnSparseClass Mat;
	Mat.load_matrix(mat_file);
	Mat.print_mat();

	Mat.diagonalize();
	Mat.save_eigen(4);	
		
// Now, the same for the Sparse version.  Find lowest 4 eigenvalues and vectors:

	SparseClass Sp;
	Sp.load_matrix(mat_file);
	Sp.print_mat();
	Sp.set_nev(4);
	
	Sp.diagonalize();
	Sp.save_eigen(4);


	return 0;
}



