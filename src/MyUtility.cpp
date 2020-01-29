/**
* @file: MyUtility.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "MyUtility.hpp"
using namespace std;
using namespace Eigen;

ullong factorial(unsigned int n){
	if (n == 0)
		return 1;
	return(n * factorial(n - 1));
}

ullong factorial(int n){
	if (n == 0)
		return 1;
	return(n * factorial(n - 1));
}

void printTiePath(const int & v, // vertex number;
	const std::vector<int> & v_tiePath, // tie-paths;
	const std::string  & fn) // output file name;
{
	int row_num = v_tiePath.size() / v;
	int ** p_Gp = new int *[row_num*v];
	for (int i = 0; i < row_num*v; ++i){
		p_Gp[i] = new int[v];
	}
	// initialize
	for (int i = 0; i < row_num*v; ++i){
		for (int j = 0; j < v; ++j){
			p_Gp[i][j] = 0;
		}
	}

	std::ofstream  of;
	of.open(fn.c_str());
	if (!of.is_open()){
		std::cout << "Error! Not open " << fn << endl;
	}
	else{
		for (int r = 0; r < row_num; ++r){
			for (int c = 0; c < v - 1; ++c){
				p_Gp[r*v + v_tiePath[r*v + c]][v_tiePath[r*v + c + 1]] = 1;
			}
		}

		for (int i = 0; i < row_num*v; ++i){
			for (int j = 0; j < v; ++j){
				std::cout << p_Gp[i][j] << ",";
				of << p_Gp[i][j] << ",";
			}
			std::cout << endl; of << endl;
		}
	}
	// release memory
	for (int i = 0; i < row_num*v; ++i){
		delete[] p_Gp[i];
	}
	delete[] p_Gp;
}

string toString(const vector<int> &p) {
	string s;

	for (vector<int>::const_iterator it = p.begin(); it != p.end(); it++) {
		s += std::to_string(*it) + " <- ";
	}

	return s;
}

void printFactorial(std::ofstream  & of, const int & n){
	ullong delta = factorial(n) / 10;
	of << "factorial(" << n << ") = " << factorial(18) << "\n"
		<< "delta = factorial(" << n << ")/10 = " << delta << "\n";
	for (int i = 0; i < 10; ++i){
		of << "[  " << i + 1 << " " << i*delta << " " << (i + 1)*delta << "  ]\n";
	}

	for (int i = 0; i < 10; ++i){
		of << i + 1 << " " << i*delta << " " << (i + 1)*delta << "\n";
	}
}


std::vector<std::vector<double> >
add(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs) {
	vector<vector<double> > result;

	result = lhs;
	for (int i = 0; i < (int)lhs.size(); i++) {
		for (int j = 0; j < (int)lhs[0].size(); j++)
			result[i][j] = lhs[i][j] + rhs[i][j];
	}

	return result;
}

std::vector<std::vector<int>> add(const std::vector<std::vector<int> >& lhs, const std::vector<std::vector<int> >& rhs){
	vector<vector<int> > result;

	result = lhs;
	for (int i = 0; i < (int)lhs.size(); i++) {
		for (int j = 0; j < (int)lhs[0].size(); j++)
			result[i][j] = lhs[i][j] + rhs[i][j];
	}
	return result;
}

std::vector<std::vector<double> >
multiplyEigen(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs){
	int dim1 = 0, dim2 = 0, dim3 = 0;

	//get dimensions of lhs and rhs
	dim1 = (int)lhs.size();
	dim2 = (int)lhs[0].size();
	dim3 = (int)rhs[0].size();
	/*
	Sure thing. You can't do the entire matrix at once, 
	because vector<vector> stores single rows in contiguous memory, 
	but successive rows may not be contiguous. 
	But you don't need to assign all elements of a row:

	std::vector<std::vector<double> > data;
	MatrixXd mat(10, 4);
	for (int i = 0; i < 10; i++)
	mat.row(i) = VectorXd::Map(&data[i][0],data[i].size);
	*/

	Eigen::Matrix2d Lh(dim1, dim2), Rh(dim2, dim3), Re(dim1, dim3);
	for (int i = 0; i < dim1; i++){
		Lh.row(i) = Eigen::VectorXd::Map(&lhs[i][0], dim2);
	}
	for (int i = 0; i < dim2; i++){
		Rh.row(i) = Eigen::VectorXd::Map(&rhs[i][0], dim3);
	}

	Re = Lh * Rh;
	vector<vector<double> > result(dim1, vector<double>(dim2, 0));
	for (int i = 0; i < dim1; i++){
		for (int j = 0; j < dim2; j++){
			result[i][j] = Re(i,j);
		}
	}
	return result;
}

std::vector<std::vector<double> >
multiply(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs) {
	int dim1, dim2, dim3;
	vector<vector<double> > result;
	vector<double> row;

	//get dimensions of lhs and rhs
	dim1 = (int)lhs.size();
	dim2 = (int)lhs[0].size();
	dim3 = (int)rhs[0].size();

	//initialize result matrix
	for (int i = 0; i < dim1; i++) {
		row = vector<double>(dim3);
		result.push_back(row);
	}

	//assign value
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim3; j++) {
			result[i][j] = 0;
			for (int k = 0; k < dim2; k++)
				result[i][j] += lhs[i][k] * rhs[k][j];
		}
	}

	return result;
}

// direct use of weights.
void unify_BT_D(std::vector<std::vector<double> >& M) {
	size_t dim1 = M.size();
	double sum;

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < i; j++) {
			sum = M[i][j] + M[j][i];
			if(sum != 0) {
				M[i][j] = M[i][j] / sum;
				M[j][i] = 1.0 - M[i][j];
			}

		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < dim1; ++i)
		M[i][i] = 0.0;
}

void unify_BT_D(Eigen::MatrixXd & M, const int & n){
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			sum = M(i, j) + M(j, i);
			if (sum != 0) {
				M(i, j) /= sum;
				M(j, i)/= sum;
			}

		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < n; ++i)
		M(i,i) = 0.0;
}

// exponential weights
void unify_BT_E(vector<vector<double> >& M) {
	size_t dim1;
	double sum;
	dim1 = M.size();

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < i; j++) {
			sum = exp(M[i][j]) + exp(M[j][i]);
			if (sum != 0) {
				M[i][j] = exp(M[i][j]) / sum;
				M[j][i] = exp(M[j][i]) / sum;
			}
		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < dim1; ++i)
		M[i][i] = 0.0;
}

/*Laplace smoothing*/
void unify_BT_L(std::vector<std::vector<double> >& M, const double & alpha){
	size_t dim1;
	double sum;
	dim1 = M.size();

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < i; j++) {
			sum = M[i][j] + M[j][i] + 2*alpha;
			if (sum != 0) {
				M[i][j] = ( alpha + M[i][j]) / sum;
				M[j][i] = ( alpha + M[j][i]) / sum;
			}
		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < dim1; ++i)
		M[i][i] = 0.0;
}

void output_matrix(std::vector<std::vector<int> > & matrix) {
	int n;

	n = (int)matrix.size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}


void
output_matrix(std::vector<std::vector<double> > &matrix) {
	int n;

	n = (int)matrix.size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}