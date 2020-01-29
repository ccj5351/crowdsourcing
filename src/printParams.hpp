/**
* @file: printParams.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_PRINT_TASK_ASSIGNMENT_PARAMETERS_H_
#define __HEADER__CROWD_SOURCING_PRINT_TASK_ASSIGNMENT_PARAMETERS_H_
#include <iostream>
#include <fstream> // std::fstream
#include <string>
#include <vector>
#include <cstdio>/*printf*/

using namespace std;
//using namespace Eigen;

void print_task_assignment_parameters(
	const bool & IsSave2File,
	ofstream & myfile,
	const int & case_idx,
	const int & n,
	const int & l,
	const vector<int> & v_m,
	const vector<int> & v_B);

#define _print_task_assignment_parameters__
#ifndef _print_task_assignment_parameters__
// example to run the function in this file;
int main(){
	//print_task_assignment_parameters

	vector<int> v_m = { 40, 60, 100, 120, 150, 180, 200 };
	vector<int> v_B = { 10, 26, 50, 100, 130, 260, 300, 350, 400 };
	int case1 = 1, n1 = 30, l1 = 2;
	int case2 = 2, n2 = 30, l2 = 3;
	int case3 = 3, n3 = 50, l3 = 2;
	int case4 = 4, n4 = 50, l4 = 3;
	bool IsSave2File = false;
	string fn = "/home/ccj/task-assignment-result.txt";
	ofstream myfile(fn);
	if (myfile.is_open()){
		print_task_assignment_parameters(IsSave2File, myfile, case1, n1, l1, v_m, v_B);
		print_task_assignment_parameters(IsSave2File, myfile, case2, n2, l2, v_m, v_B);
		print_task_assignment_parameters(IsSave2File, myfile, case3, n3, l3, v_m, v_B);
		print_task_assignment_parameters(IsSave2File, myfile, case4, n4, l4, v_m, v_B);
	}
	myfile.close();

	return 0;
}
#endif
#endif