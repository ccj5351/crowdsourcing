/**
* @file: transition.cpp
* @brief: to compute the transitivity closure.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "transition.hpp"
#include "MyUtility.hpp"
#include "Experiment.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
using namespace Eigen;



std::vector<std::vector<double> >
generateTransitionMatrix(const std::string & path, const int &n) {
	vector<vector<double> > M, result;
	vector<vector<vector<double> > > transitions;
	ifstream ifs;
	istringstream iss;
	double count;
	string line;

	//read matrix M from file
	ifs.open(path);
	for (int i = 0; i < n; i++) {
		M.push_back(vector<double>(n, .0));
		getline(ifs, line);
		iss.str(line);
		for (int j = 0; j < n; j++) {
			iss >> count;
			M[i][j] = count;
		}
	}

	unify_BT_D(M);

	transitions.push_back(M);
	//for i = 2 to n-1
	for (int i = 2; i < n; i++) {
		//do multiplication and unify
		result = multiply(M, transitions[i - 2]);
		// unify_BT_D(result);
		//store
		transitions.push_back(result);
	}

	//do sum and unify
	result = transitions[0];
	for (int i = 1; (int)i < transitions.size(); i++) {
		result = add(result, transitions[i]);
	}


	/*
	for (int i=0; i<n; i++) {
	for (int j=0; j<n; j++)
	cout << result[i][j] << " ";
	cout << endl;
	}
	*/

	return result;
}


std::vector<std::vector<double> >
generateTransitionMatrix(const std::vector<std::vector<double> >& matrix){

	vector<vector<double> > M, result;
	vector<vector<vector<double> > > transitions;
	int n;

	M = matrix;

	n = (int)M.size();

	unify_BT_D(M);

	transitions.push_back(M);
	//for i = 2 to n-1
	for (int i = 2; i < n; i++) {
		//do multiplication and unify
		result = multiply(transitions[i - 2], M);
		// Pay attention : matrix a*b != b*a;
		// the order matter for matrix multiplication.
		//result = multiply(M, transitions[i - 2]);
		unify_BT_D(result);
		/*cout << "i = " << i << endl;
		print_2D_vector<double>(n, result);*/
		
		//store
		transitions.push_back(result);
	}

	//do sum and unify
	result = transitions[0];
	for (int i = 1; (int)i < transitions.size(); i++) {
		result = add(result, transitions[i]);
	}

	unify_BT_D(result);
	return result;

}

// Due to the O(n^3) complexity of matrix multiplication,
// here we use the Eigen library to do the matrix multiplication to obtain
// efficient computation, even though it will consume CPU too much.
std::vector<std::vector<double> >
generateTransitionMatrix_Eigen(const std::vector<std::vector<double>> & matrix){

	int n = (int)matrix.size();
	vector<vector<double>> M(matrix);
	unify_BT_D(M);

	MatrixXd cur(n, n), m(n, n), sum_(n, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j){
			cur(i, j) = M[i][j];
			m(i, j) = M[i][j];
			sum_(i, j) = M[i][j];
		}
	for (int i = 0; i < n - 2; i++) {
		cur = cur * m;
		// Option 1:
		unify_BT_D(cur, n);
		/*
		cout << "i = " << i << endl;
		print_2D_Eigen(n, cur);
		*/

		sum_ += cur;
		// Option 2:
		//unify_BT_D(cur);
		/*cout << " i = " << i << " , cur (2, 5) = " << cur(2, 5) << ", 
		  sum (2, 5) = " << sum_(2, 5) 
			<<",  /= " 
			<< sum_(2, 5) / (sum_(2, 5) + sum_(5,2)) << endl;*/
	}
	//do sum and unify
	unify_BT_D(sum_, n);

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			M[i][j] = sum_(i, j);
		M[i][i] = 0;
	}
	return M;
}