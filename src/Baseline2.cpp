/**
* @file: Baseline2.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include "Baseline2.hpp"
using namespace std;

//this is the implementation of Condorcet Fusion with QuickSort aggregation
//the input should be voting matrix with integers as elements.
/*see the paper "Condorcet Fusion for Improved Retrieval"
 *2002 year, CIKM.
 */


// transitive closure matrix;
int* MajorityVoting(const int & vertexNum, const vector<vector<double>> & G_TC){
 
    //calculate the weight of each node
	vector<std::pair<int, double>> v_weight(vertexNum);
    for ( int i = 0 ; i < vertexNum ; ++ i){
		v_weight[i] = make_pair(i, .0);
    }
    
	// in ----> +;
	// out ---> -;
    for ( int i = 0 ; i < vertexNum ; ++ i){
        for ( int j = 0 ; j < vertexNum ; ++ j){
            v_weight[i].second += G_TC[j][i];
        }
    }

    for ( int i = 0 ; i < vertexNum ; ++ i){
        for ( int j = 0 ; j < vertexNum ; ++ j){
            v_weight[i].second -= G_TC[i][j];
        }
    }
    
	std::sort(v_weight.begin(), v_weight.end(), 
		[](const std::pair<int, double> & arg1,
		const std::pair<int, double> & arg2){return arg1.second <= arg2.second; });
	int * p_rank = new int[vertexNum];
    for ( int i = 0 ; i < vertexNum ; ++ i){
		p_rank[i] = v_weight[i].first;
    }

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	cout << "Baseline ranking result is :\n";
    for ( int i = 0 ; i < vertexNum ; ++ i){
        cout << p_rank[i] << ",  " ;
    }
	cout << std::endl;
#endif

	vector<std::pair<int, double>>().swap(v_weight);
    return p_rank ;
}


// input matrix: voting result;
int* MajorityVoting(const int & vertexNum, const vector<vector<int>> & Voting){

	//calculate the weight of each node
	vector<std::pair<int, int>> v_scores(vertexNum);
	for (int i = 0; i < vertexNum; ++i){
		v_scores[i] = make_pair(i, 0);
	}

	// in ----> +;
	// out ---> -;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_scores[i].second += Voting[j][i];
		}
	}

	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_scores[i].second -= Voting[i][j];
		}
	}

	//sort works with lambda function in c++0x/c++11;
	std::sort(v_scores.begin(), v_scores.end(), 
		[](const std::pair<int, int> & arg1,
		const std::pair<int, int> & arg2){return arg1.second <= arg2.second; });
	int * p_rank = new int[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		p_rank[i] = v_scores[i].first;
	}

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	cout << "Baseline ranking result is :\n";
	for (int i = 0; i < vertexNum; ++i){
		cout << p_rank[i] << ",  ";
	}
	cout << std::endl;
#endif

	vector<std::pair<int, int>>().swap(v_scores);
	return p_rank;
}

/*high rank means better priority*/
bool IsPrefer_J_Than_I(const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M, // each voter's voting result;
	// that is: if voter prefers i than j, then edge (i, j) exists, i.e., v[i][j] = 1, v[j][i] = 0;
	const int & i, // object i
	const int & j // object j
	){
	int count = 0;
	for (int k = 0; k < w; ++k){
		if (v_M[k][i][j] == 1)
			count--;
		if (v_M[k][j][i] == 1)
			count++;
	}
	/*if count > 0, prefers j than i*/
	if (count >= 0)
		return true;
	else
		return false;
}


int partition(int *arr, const int left, const int right,
	const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	) {
	const int mid = left + (right - left) / 2;
	const int pivot = arr[mid];
	// move the mid point value to the front.
	std::swap(arr[mid], arr[left]);
	int i = left + 1;
	int j = right;
	while (i <= j) {
		while (i <= j && IsPrefer_J_Than_I(w, v_M, arr[i], pivot)){
			i++;
		}

		while (i <= j && IsPrefer_J_Than_I(w, v_M, pivot, arr[j])) {
			j--;
		}

		if (i < j) {
			std::swap(arr[i], arr[j]);
		}
	}
	std::swap(arr[i - 1], arr[left]);
	return i - 1;
}

/*the final sorted result will be saved in the "arr"*/
void CondorcetFuseSort(int *arr, const int left, const int right, const int sz,
	const int & w, // each of the w search systems, or w voters; 
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	){
	if (left >= right) {
		return;
	}
	int part = partition(arr, left, right, w, v_M);
	//std::cout << "QSC:" << left << "," << right << " part=" << part << "\n";
	//print(arr, sz);
	CondorcetFuseSort(arr, left, part - 1, sz, w, v_M);
	CondorcetFuseSort(arr, part + 1, right, sz, w, v_M);
}
