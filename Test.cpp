#include "Test.h"

void test_load_UnSparseClass(){
	
	string mat_file = "matrix.dat";

	UnSparseClass Mat;
	Mat.load_matrix(mat_file);
	Mat.print_mat();

	cout << "\n Copy UnSparseClass \n";
	UnSparseClass *Mat_copy = new UnSparseClass(Mat);
	Mat_copy->print_mat();

	return;
}

void test_load_SparseClass(){

	string mat_file = "matrix.dat";

	SparseClass Mat;
	Mat.load_matrix(mat_file);
	Mat.print_mat;
	Mat2.show_compression();

	cout << "\n Copy SparseClass \n";
	SparseClass *Mat_copy = new SparseClass(Mat);
	Mat_copy-->print_mat();

	return;
}

