/**
* @file: Baseline2.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_BASE_LINE2_H_
#define __HEADER__CROWD_SOURCING_BASE_LINE2_H_

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>

template <typename T>
bool ascendSorting(const std::pair<int, T> & arg1,
	const std::pair<int, T> & arg2){
	return (arg1.second <= arg2.second);
}

// transitive closure matrix;
int * MajorityVoting(const int & vertexNum, 
	const std::vector<std::vector<double>> & G_TC);
int * MajorityVoting(const int & vertexNum, 
	const std::vector<std::vector<int>> & Voting);

/*high rank means better priority*/
bool IsPrefer_J_Than_I(const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M, // each voter's voting result;
	// that is: if voter prefers i than j, then edge (i, j) exists, i.e., v[i][j] = 1, v[j][i] = 0;
	const int & i, // object i
	const int & j // object j
	);

int partition(int *arr, const int left, const int right,
	const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	);

/*the final sorted result will be saved in the "arr"*/
void CondorcetFuseSort(int *arr, const int left, const int right, const int sz,
	const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	);

#endif