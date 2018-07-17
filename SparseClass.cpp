#include "SparseClass.h"

extern "C" {
	void dsaupd_( int*, char*, int*, char*, int*, double*, double*, int*, 
				  double*, int*, int*, int*, double*, double*, int*, int*);

	void dseupd_( int*, char*, int*, double*, double*, int*, double*, char*, 
				  int*, char*, int*, double*, double*, int*, double*,
				  int*, int*, int*, double*, double*, int*, int*);
}

SparseClass::SparseClass(const SparseClass &sp){		// Copy Constructor	
	
	dim = sp.dim;
	num_non_zero = sp.num_non_zero;
	zero_tol = sp.zero_tol;

	int i;
	diagonal = new double [dim];
    for(i = 0; i< dim; i++)
		diagonal[i] = sp.diagonal[i];

	first_in_row = new int [dim+1];
	for(i=0; i< dim+1; i++)
		first_in_row[i] = sp.first_in_row[i];

	values = sp.values;
	column_index = sp.column_index;
	allocate_evals();
	
	return;
}

SparseClass &SparseClass:: operator=(const SparseClass &sp){ 	// Copy assignment constructor
	
	if( this != &sp){
		// deallocate memory
		delete [] diagonal;
		delete [] first_in_row;
		values.clear();
		column_index.clear();
	
		// Now, store new variables	
		dim = sp.dim;
		num_non_zero = sp.num_non_zero;
		zero_tol = sp.zero_tol;

		int i;
		diagonal = new double [dim];
	    for(i = 0; i< dim; i++)
			diagonal[i] = sp.diagonal[i];

		first_in_row = new int [dim+1];
		for(i=0; i< dim+1; i++)
			first_in_row[i] = sp.first_in_row[i];

		values = sp.values;
		column_index = sp.column_index;
		allocate_evals();
	}
	return (*this);
}	

SparseClass::~SparseClass(){ 	// Destructor
	delete [] diagonal;
	delete [] first_in_row;
	values.clear();
	column_index.clear();

	delete [] eigenvalues;
	for(int i = 0; i < dim; i++)
		delete [] eigenvectors[i];
	delete [] eigenvectors;

	return;
}

void SparseClass::load_matrix(string infile){  // virtual function
	
	filename = infile;
	find_dim();
	cout << "/nSparse load_matrix()\n";

	// By default, set num_ev = dim-1; 
	// it can be changed with set_nev() 
	num_ev = dim-1;

	// First, read matrix in from file into temporary storage:
	
    ifstream file(filename);
	int i, j;
	double ** temp = new double*[dim];

	for(i=0; i<dim; i++){
		temp[i] = new double [dim];
		for(j = 0; j < dim; j++){
			file >> temp[i][j];
		}	
	}
	file.close();

	cout << "Matrix in temp storage\n";
	// 1. Load diagonal
    diagonal = new double [dim];
	for(i = 0; i< dim; i++)
		diagonal[i] = temp[i][i];

	cout << "diagonal stored\n";
	// 2. Load values, first_in_row, column_index.  As a reminder:
	//    vector<double> values[i] ==  i'th non-zero matrix element of lower triangle, row-major order
	//    int * first_in_row[k]  ==  index in values of first element of row k, k = 0, 1, ..., dim
	//    vector<int> column_index[i] == column index of values[i] in matrix 
	
	// 
	first_in_row = new int[dim+1]; 
	first_in_row[0] =  0;   	// Row 0 has no elements before the diagonal;
	bool it_is_first = true;   	// True until first non-zero in row is found, then
          						// remains false until end of row; then resets
								
	int last_in_values = 0;		// within row i, holds index in values of most
								// recently added element.
								// at end of a row, last_in_values+1 will be
								// added to the next first_in_row slot.
	

	for( i=1; i<dim; i++){
		for( j = 0; j < i; j++){
			if( temp[i][j] > zero_tol || temp[i][j] < zero_tol){ // if temp[i][j] != 0
				values.push_back(temp[i][j]);
				column_index.push_back(j);
				last_in_values++;
				if(it_is_first){
					first_in_row[i] = values.size()-1;
					it_is_first = false;
				}	
			}	
		}
		if(it_is_first)  // if there is no non-zero element before the diagonal...
			first_in_row[i] = last_in_values;
		it_is_first = true;  // resets at end of row
	}

	num_non_zero = values.size();	
	first_in_row[dim] = num_non_zero;

	allocate_evals();

	cout << "delete temp[][]\n";
	for(i = 0; i<dim; i++)
		delete [] temp[i];
	delete [] temp;
	
	cout << "end load_matrix()\n";
	return;
}	

void SparseClass::show_compression(){
	cout << "In sparse compression, the matrix is described by 4 vectors:\n";
	cout << "dimension = " << dim << endl;
	cout << "num_non_zeros = " << num_non_zero << endl;

	int i;
	cout << "\n diagonal[] = ";
	for(i = 0; i<dim; i++)
		cout << " " << diagonal[i] ;
	cout << "\n values[] = ";
	for(i = 0; i<num_non_zero; i++)
		cout << " " << values[i] ;
	cout << "\n first_in_row[] = ";
	for(i = 0; i<dim+1; i++)
		cout << " " << first_in_row[i] ;
	cout << "\n column_index[] = ";
	for(i = 0; i<num_non_zero; i++)
		cout << " " << column_index[i] ;
	
}


void SparseClass::print_mat(){
	int i, j;
	cout << "upper left 4x4 block of matrix: \n";
	for(i=0; i<4; i++){
		for(j=0; j<4; j++)
			cout << " " << get_elem_ij(i, j);
		cout << endl;
	}
	return;
}


// Return row i, column j of matrix stored in compressed format
double SparseClass::get_elem_ij(int i, int j){   // virtual function
	double element = 0;

	int first=0; // index in values of first element in row i
	int last=0;  // index in values[] of last element in row i

	int k;  // look at column_index[k] between k = first and k=last to see 
			// when and if column_index[k] == j
			
	if(i==j)
		element = diagonal[i];

	else if(i > j){  // Lower left triangle, row > column.
		first = first_in_row[i];
		last = first_in_row[i+1] - 1;  

		if( (last-first) < 0)  	// There are no elements in that row, all are zero
			element = 0;
		else{	 				// There is at least one non-zero element in the row
			element = 0;
			for( k = first; k <= last; k++){
				if(column_index[k] == j)
					element = values[k];
			}
		}
	}
	else if(i < j){  // Upper right triangle, row < column
		// Same as before, but now look at row = j, column = i.
		element = get_elem_ij(j, i);
	}
	
	return element;
}

void SparseClass::set_nev(int num_evals){
	num_ev = min(num_evals, dim-1);
	return;
}

void SparseClass::diagonalize(){
	
	int n = dim;
	dsaupd(n, num_ev, eigenvalues, eigenvectors);
	return;
}

/* ------------------------------------------------------------------
 * ARPACK interface functions:  the following two functions,
 * dsaupd and dseupd, are modified from a version
 * published online by Scot Shaw (30 August, 1999)
 * 
 * For these functions to work, there must be a function defined as followe:
 * av( int n, double *in, double *out)
 * which defines how the user-defined matrix acts on a generic vector (*in) to 
 * produce a resulting vector (*out)
 *
 * This function is known as the reverse-communication interface in ARPACK.
 */
void SparseClass::dsaupd(int n, int nev, double *Evals)
{
  int ido = 0; /* Initialization of the reverse communication
                parameter. */
  
  char bmat[2] = "I"; /* Specifies that the right hand side matrix
                       should be the identity matrix; this makes
                       the problem a standard eigenvalue problem.
                       Setting bmat = "G" would have us solve the
                       problem Av = lBv (this would involve using
                       some other programs from BLAS, however). */
  
  char which[3] = "SM"; /* Ask for the nev eigenvalues of smallest
                         magnitude.  The possible options are
                         LM: largest magnitude
                         SM: smallest magnitude
                         LA: largest real component
                         SA: smallest real compoent
                         LI: largest imaginary component
                         SI: smallest imaginary component */
  
  double tol = 0.0; /* Sets the tolerance; tol<=0 specifies
                     machine precision */
  
  double *resid;
  resid = new double[n];
  
  int ncv = 4*nev; /* The largest number of basis vectors that will
                    be used in the Implicitly Restarted Arnoldi
                    Process.  Work per major iteration is
                    proportional to N*NCV*NCV. */
  if (ncv>n) ncv = n;
  
  double *v;
  int ldv = n;
  v = new double[ldv*ncv];
  
  int *iparam;
  iparam = new int[11]; /* An array used to pass information to the routines
                         about their functional modes. */
  iparam[0] = 1;   // Specifies the shift strategy (1->exact)
  iparam[2] = 3*n; // Maximum number of iterations
  iparam[6] = 1;   /* Sets the mode of dsaupd.
                    1 is exact shifting,
                    2 is user-supplied shifts,
                    3 is shift-invert mode,
                    4 is buckling mode,
                    5 is Cayley mode. */
  
  int *ipntr;
  ipntr = new int[11]; /* Indicates the locations in the work array workd
                        where the input and output vectors in the
                        callback routine are located. */
  
  double *workd;
  workd = new double[3*n];
  
  double *workl;
  workl = new double[ncv*(ncv+8)];
  
  int lworkl = ncv*(ncv+8); /* Length of the workl array */
  
  int info = 0; /* Passes convergence information out of the iteration
                 routine. */
  
  int rvec = 0; /* Specifies that eigenvectors should not be calculated */
  
  int *select;
  select = new int[ncv];
  double *d;
  d = new double[2*ncv]; /* This vector will return the eigenvalues from
                          the second routine, dseupd. */
  double sigma;
  int ierr;
  char All[3] = {'A','l','l'};
  
  /* Here we enter the main loop where the calculations are
   performed.  The communication parameter ido tells us when
   the desired tolerance is reached, and at that point we exit
   and extract the solutions. */
  
  do {
    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, &info);
    
    if ((ido==1)||(ido==-1)) av(n, workd+ipntr[0]-1, workd+ipntr[1]-1);
  } while ((ido==1)||(ido==-1));
  
  /* From those results, the eigenvalues and vectors are
   extracted. */
  
  if (info<0) {
    cout << "Error with dsaupd, info = " << info << "\n";
    cout << "Check documentation in dsaupd\n\n";
  } else {
    dseupd_(&rvec, All, select, d, v, &ldv, &sigma, bmat,
            &n, which, &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &ierr);
    
    if (ierr!=0) {
      cout << "Error with dseupd, info = " << ierr << "\n";
      cout << "Check the documentation of dseupd.\n\n";
    } else if (info==1) {
      cout << "Maximum number of iterations reached.\n\n";
    } else if (info==3) {
      cout << "No shifts could be applied during implicit\n";
      cout << "Arnoldi update, try increasing NCV.\n\n";
    }
    
    /* Before exiting, we copy the solution information over to
     the arrays of the calling program, then clean up the
     memory used by this routine.  For some reason, when I
     don't find the eigenvectors I need to reverse the order of
     the values. */
    
    int i;
    for (i=0; i<nev; i++) Evals[i] = d[nev-1-i];
    
    delete resid;
    delete v;
    delete iparam;
    delete ipntr;
    delete workd;
    delete workl;
    delete select;
    delete d;
  }
  return;
}


void SparseClass::dsaupd(int n, int nev, double *Evals, double **Evecs)
{
  int ido = 0;
  char bmat[2] = "I";
  char which[3] = "SM";
  double tol = 0.0;
  double *resid;
  resid = new double[n];
  int ncv = 4*nev;	// info = -3: NCV must be greater than NEV and less than or equal to N
  if (ncv>n) ncv = n;
  double *v;
  int ldv = n;
  v = new double[ldv*ncv];
  int *iparam;
  iparam = new int[11];
  iparam[0] = 1;
  iparam[2] = 3*n;
  iparam[6] = 1;
  int *ipntr;
  ipntr = new int[11];
  double *workd;
  workd = new double[3*n];
  double *workl;
  workl = new double[ncv*(ncv+8)];
  int lworkl = ncv*(ncv+8);
  int info = 0;
  int rvec = 1;  // Changed from above
  int *select;
  select = new int[ncv];
  double *d;
  d = new double[2*ncv];
  double sigma;
  int ierr;
  char All[3] = {'A','l','l'};
  
  do {
    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, &info);
    
    if ((ido==1)||(ido==-1)) av(n, workd+ipntr[0]-1, workd+ipntr[1]-1);
  } while ((ido==1)||(ido==-1));
  
  if (info<0) {
    cout << "Error with dsaupd, info = " << info << "\n";
    cout << "Check documentation in dsaupd\n\n";
  } else {
    dseupd_(&rvec, All, select, d, v, &ldv, &sigma, bmat,
            &n, which, &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &ierr);
    
    if (ierr!=0) {
      cout << "Error with dseupd, info = " << ierr << "\n";
      cout << "Check the documentation of dseupd.\n\n";
    } else if (info==1) {
      cout << "Maximum number of iterations reached.\n\n";
    } else if (info==3) {
      cout << "No shifts could be applied during implicit\n";
      cout << "Arnoldi update, try increasing NCV.\n\n";
    }
    
    int i, j;
    for (i=0; i<nev; i++) Evals[i] = d[i];
    for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[j][i] = v[i*n+j];
    
    delete resid;
    delete v;
    delete iparam;
    delete ipntr;
    delete workd;
    delete workl;
    delete select;
    delete d;
  }
  return;
}

/*
 *	Reverse communication interface av().
 *
 *	This function defines the product of multiplying our user-defined sparse
 *	matrix times a column vector, 
 *	array[][] * in[] = out[]
 *	In this example, the matrix has been stored
 *	in memory using a packed storage scheme.  However, for many projects, it is
 *	more efficient to calculate the matrix elements "on-the-fly" for this 
 *	Matrix-on-vector multiplication function,  instead of storing them in memory
 *	to be looked up.
 *
 *  
 */ 

void SparseClass::av(int n, double *in, double *out){
	int i;  

	// First, add the diagonal elements: A[i][i] * in[i] = out[i]	
	for( i = 0; i < n; i++)
		out[i] = diagonal[i]*in[i];

	// Now, move through values[] and add lower left and upper right triangle
	// parts to appropriate out[j]
	// Lower left triangle: out[j] = M[i][j]*in[j]
	// UpperRight triangle: out[j] = M[j][i]*in[j]
	//      ...using Einstein summation convention: sum over repeated indices (j)

	// First, find where values lie in matrix: find row, column for that element: 
	int j = 0;  // index in first_in_row[], start at zero
	int k = 0;  // index in Values, 0 <= k < num_nun_zero

	/* The next while loops will increment k and j, 
	 * When first_in_row[j] <= k < first_in_row[j+1], we're in row j of the matrix
	 * increment j until end of row
	 * add lower triangle and upper triangle product terms to out[], then 
	 * increment k to move to next non-zero matrix element
	 */
	
	while( k < num_non_zero){
		while( k >= first_in_row[j+1]){
			j++;
		}	// now, j = row of values[k], column = column_index[k]
		out[j] += values[k]*in[column_index[k]];
		out[column_index[k]] += values[k] * in[j];
		k++;
	}

	return;
}



void SparseClass::save_eigen(int how_many){
	
	ofstream printfile;
	printfile.open("evals_Arpack.dat");
	int i;
	for( i = 0;  i < min(num_ev, how_many); i++){
		printfile << eigenvalues[i] << endl;
	}
	printfile.close();

	printfile.open("evecs_Arpack.dat");
	int j;
	for(i = 0; i < min(num_ev, how_many); i++){
		for(j = 0; j < dim; j++){
			printfile << eigenvectors[i][j] << " ";
		}
		printfile << endl;
	}
	printfile.close();
	return;
}


