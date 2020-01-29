/**
* @file: task_assign.hpp
* @brief: task graph generation.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_TASK_GENERATION_H_
#define __HEADER__CROWD_SOURCING_TASK_GENERATION_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <map>

typedef std::pair<int, int> EDGE;
using namespace std;

void deletElement(vector<int> &SetV, int element);
bool deletEdge(vector<EDGE> &SetV, EDGE element);
void printMatrix(const int & vertexNum, bool ** matrix_Gt);

bool isDegreeSuccess(const int & vertexNum, bool ** matrix_Gt, const int & deg);

bool assign_task_graph(
	const int & n, // vertex num
	const int & d, // degree
	bool ** g_t
	);

bool assign_task_graph_down_2_d(
	const int & n, // vertex num
	const int & d, // degree
	bool ** g_t,
	ofstream & of
	);

#endif