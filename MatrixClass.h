#ifndef MATRIXCLASS_H
#define MATRIXCLASS_H

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class MatrixClass{

	public:
		MatrixClass(){};

		virtual ~MatrixClass()= default;	// Abstract class: MatrixClass on its own isn't useful

		virtual void print_mat(){ cout << "MatrixClass print_mat() does nothing\n";}

		void allocate_evals();
		virtual void load_matrix(string infile) = 0;              
		virtual double get_elem_ij(int i, int j) = 0; 

		virtual void diagonalize(){};	
		virtual void save_eigen(int how_many){cout << how_many << endl;}


	protected:
		string filename;	
		int dim;
		void find_dim();  														// find_dim() sets dimension of matrix

		double *eigenvalues;
		double **eigenvectors;


	private:
		
};

#endif
