/**
* @file: simulation.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_SIMULATION_H_
#define __HEADER__CROWD_SOURCING_SIMULATION_H_
#include<iostream>
#include<windows.h>
#include <stdlib.h>
#include <time.h>  
#include <fstream>
#include <string>
#include<cstring>
using namespace std;
//read the edge information matrix
bool scanMatrix(double **  matrix, const int & n,
	const std::string & filename);
bool scanMatrix2(double **  matrix, const int & n,
	const std::string & filename);

//generate uniform number
void generate_uniform_num(int n, double* num);

//uniform distribution simulation
void simulate_order(double ** task, double ** order_matrix, 
	double *uniform_num, int n, int edge_num);
void clear_matrix(double ** matrix, int n);
void getSimuData(const int & edgeNum, // e.g. edge_num = 20;
	const int & vertexNum, // n = 10;
	const std::string & inputFileName,
	//e.g., "C:\\Users\\Mobrick\\Desktop\\simulation\\matrix.txt"
	// matrix rows are separated by ";" for MATLAB readable-data;
	// matrix rows are separated by "," for C++ readable-data;
	const std::string & outputFileM, // for store MATLAB data
	// all the data in the matrix is output 
	// as separated rows.
	const std::string & outputFileC, // for store C++ data;
	// all the data in the matrix is output 
	// as a whole line.
	const std::string & outputFileC1 // for store C++ data in one line;
	);
#endif