#include "UnSparseClass.h"

extern "C" {
	extern int dspev_( char*, char*, int*, double*, double*, 
		    double*, int*, double*, int*);
		   	
}

UnSparseClass::UnSparseClass(const UnSparseClass &un){	// Copy constructor

	dim = un.dim;
	array = new double [dim*(dim+1)/2];
	for(int i=0; i< (dim*(dim+1))/2; i++){
		array[i] = un.array[i];
	}

	allocate_evals();
	return;
}

UnSparseClass &UnSparseClass:: operator=(const UnSparseClass &un){  // Copy assignment constructor
	if( this != &un){
		// deallocate memory
		int i;
		delete [] array;

		// store new variables
		dim = un.dim;
		array = new double [dim*(dim+1)/2];
		for( i=0; i< (dim*(dim+1))/2; i++){
			array[i] = un.array[i];
		}
		
		allocate_evals();
    }
	return (*this);
}

UnSparseClass::~UnSparseClass(){	// Destructor
	delete [] array;
	delete [] eigenvalues;
	
	for( int i = 0; i<dim; i++)
			delete [] eigenvectors[i];
	delete [] eigenvectors;

	return;
}

void UnSparseClass::load_matrix(string infile){
   
	int i,j;

	filename = infile;
	find_dim();

	ifstream file(filename);
	double **temp  = new double*[dim];
	for(i=0; i<dim; i++){
		temp[i] = new double [dim];
        for(j=0; j<dim; j++)
			file >> temp[i][j];	
	}
	file.close();

	array = new double[(dim*(dim+1))/2];
  	// aij stored in array[ i + j(j+1))/2, i <= j, where 
	//  Fortran uses Column ordered storage, a11, a12, a22, a13, a23, a33	
	for ( j = 0; j < dim; j++){
		for( i = 0; i <= j; i++){
			array[ i + (j*(j+1))/2 ] = temp[i][j];
		}
	}
	allocate_evals();

	for( i = 0; i < dim; i++)
		delete [] temp[i];
	delete [] temp;
	return;
}

void UnSparseClass::print_mat(){
	int i, j;
	cout << "upper left 4x4 block of matrix: \n";

	for(i=0; i<4; i++){
		for(j=0; j<4; j++)
			cout << " " << get_elem_ij(i,j);
		cout << endl;
	}
	return;
}

double UnSparseClass::get_elem_ij(int i, int j){
		
	if( i <= j)
		return	array[ i + (j*(j+1))/2 ] ;
	else
		return get_elem_ij(j, i);
}

// The original array in upper packed storage will be destroyed in this
// procedure
void UnSparseClass::diagonalize(){
	

	cout << "eigenvalues[0] = " << eigenvalues[0] << endl;
	

	char jobz = 'V';  // Compute eigenvalues
	char uplo = 'U';  // Upper packed storage
	int n = dim;

	int ldz = dim;
	double *work = new double[3*dim];
	int info;

	cout << "diagonalize()\n";

	// In Fortran, all variables are passed by reference into subroutines.
	// *C++ requires an underscore following Fortran function calls.
	//
	dspev_( &jobz, &uplo, &n, &array[0], &eigenvalues[0], 
		    &eigenvectors[0][0], &ldz, &work[0], &info);

	cout << "dspev_ done\n";
	cout << "info = " << info << endl;
	cout << "eigenvalues[0] = " << eigenvalues[0];
	return;

}


void UnSparseClass::save_eigen(int how_many){
	
	ofstream printfile;
	printfile.open("evals_Lapack.dat");
	int i;
	for( i = 0; i < min(dim, how_many); i++){
		printfile << eigenvalues[i] << endl;		
	}
	printfile.close();

	printfile.open("evecs_Lapack.dat"); 
	// i'th row in output file is the i'th eigenvector 
	// (double-checked row/column ordering by comparing to 
	// Mathematica output)
	int j;
	for( i = 0; i< min(dim, how_many); i++){
		for( j = 0; j < dim; j++){
			printfile << eigenvectors[i][j] << "  ";
		}
		printfile << endl;
	}
	printfile.close();
	return;
}



