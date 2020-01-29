/**
* @file: MyUtility.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_PROJECT_UTILITY_H_
#define __HEADER__CROWD_SOURCING_PROJECT_UTILITY_H_
#include <iostream>
#include <fstream> // std::fstream
#include <vector>
#include <string>
#include <map> /*ordered map*/
#include <unordered_map> /*unordered map*/
#include <set>
#include <algorithm> /*Permutation*/
#include <functional>
#include <array> /* array */
#include <cmath> /*pow etc*/
#include <cstdio>/*printf*/
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <limits.h>
#include <cstring> /* memset */
#include <queue>
#include <set>
#include <vector>
#include <random>  /*std::uniform_int_distribution*/
#include <Eigen/Dense>
//#include "CrowdBT.h" // Crowd_BT
#include <boost/math/special_functions/digamma.hpp>

typedef unsigned long ulong;
typedef unsigned long long ullong;
ullong factorial(unsigned int n);
ullong factorial(int n);
void printTiePath(const int & v, // vertex number;
	const std::vector<int> & v_tiePath, // tie-paths;
	const std::string  & fn); // output file name;

std::string toString(const std::vector<int> &p);

void printFactorial(std::ofstream  & of, const int & n);


std::vector<std::vector<double> > add(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs);

std::vector<std::vector<int>> add(const std::vector<std::vector<int> >& lhs, const std::vector<std::vector<int> >& rhs);

std::vector<std::vector<double> > multiply(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs);
std::vector<std::vector<double> > multiplyEigen(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs);

void unify_BT_D(Eigen::MatrixXd & M, const int & n);/*direct use of voting results*/
void unify_BT_D(std::vector<std::vector<double> >& M);/*direct use of voting results*/
void unify_BT_E(std::vector<std::vector<double> >& M);/*using the exponentials of voting results*/
void unify_BT_L(std::vector<std::vector<double> >& M, const double & alpha = 1.0);/*Laplace smoothing*/


void
output_matrix(std::vector<std::vector<int> > &matrix);

void
output_matrix(std::vector<std::vector<double> > &matrix);

template<typename T>
void output_matrix(std::vector<std::vector<T> > &matrix, std::ofstream & of1){
	for (int i = 0; i < matrix.size(); ++i){
		for (int j = 0; j < matrix[i].size(); ++j){
			of1 << matrix[i][j];
			if (matrix[i].size() - 1 == j)
				of1 << "\n";
			else of1 << ",";
		}
	}
}

// print boundary;
template<typename T>
void printBoundary(const T & val, const int & vertexNum, std::ofstream & of_csv){
	for (auto i = 0; i < vertexNum; ++i)
		of_csv << val << ",";
	of_csv << endl;
}

template <typename T, typename T2>
bool isHP(T ** tc, T2 * v, const int & size){

		for (int i = 0; i < size - 1; i++)
		{
			if (tc[v[i]][v[i + 1]] == 0){
				/*std::cout << " (" << v[i] << ", " << v[i + 1] <<
					") = 0, ";*/
				return false;
			}
		}
		return true;
}




#endif
