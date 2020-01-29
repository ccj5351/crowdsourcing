/**
* @file: heurisHP.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_HEURISTIC_HP_H_
#define __HEADER__CROWD_SOURCING_HEURISTIC_HP_H_

#include "path.hpp"
#include <iostream>
#include <vector>
#include <cstdio> /*printf()*/
#include <map>
#include <string>
#include <algorithm>

using namespace std;
/* This algorithm is proposed in this paper.
* A HEURISTIC METHOD FOR THE DETERMINATION
* OF A HAMILTONIAN CIRCUIT IN A GRAPH.
* see http://journals.cambridge.org/download.php?file=%2FANZ%2FANZ28_03%2FS0334270000005439a.pdf&code=2fefdbb07e26bccf0e7fa9d615b280fd
* Note:
* Since we want to find an HP instead of HC,
* thus our algorithm implementation does not
* involve the so called anchor of a matrix.
*/

class NormalFormMatrix4HP{
private:
#define Large_Value 200 /*used to initialize min_ value*/
	enum idx_state : int  { S0 = 0, S1 = 1, S2 = 2, Error = 3};
	double ** g_p;
	const int verNum;
	
	// vector element idx is heading idx,
	// the corresponding value is the vertex idx.
	vector<int> v_heading;

	// vector element idx is node idx,
	// the corresponding value is the heading idx.
	vector<int> v_node;

	// normal path made up of the off-diagonal elements
	// of a matrix A = [a(i,j)], i.e., made up of (n-1) 
	// a(i,i+1) elements such that 0 <= i <= n-2.
	vector<int> v_normal_path;

	// std::map[], element-access complexity is O(logn).
	// therefore, we use std::vector here, since its
	// element-access complexity is O(1).
	// pair::first -- outgoing neighbors of node v;
	// pair::second -- disturb degree of outgoing neighbors of node v;
	vector<vector<std::pair<int, int>>> v_out_neibor_distur;
	
public:
	NormalFormMatrix4HP(double ** g,
		const int & v);
	bool isNormalFormMat(int & c);
	// interchanging column i and column j,
	// and the corresponding row i and row j.
	void swap_rows_clos(const int & i, const int & j);
	
	int least_disturbed_node(const int & h_i);

	// return how many non-zero nodes have been removed
	// after interchanging row/column idx 
	// and the candidate candiIdx.
	int DisturbDegree(
		const int & h_i, /*the column/row heading h_i,
						 corresponding to the current zero element
						 encountered on the normal path.*/
		const int & h_j); /*candidate heading h_j*/
	// return the node index, given a heading index.
	inline int getNodeIdx(const int & h_i)
	{
		return v_heading[h_i];
	}
	// return the heading index, given a node index.
	inline int getHeadingIdx(const int & v_i)
	{
		return v_node[v_i];
	}
	void showMatrix();
	void showMatrix(const double ** g_p_);
	inline idx_state getIdxState(const int & h){
		if (0 == h)
			return S0;
		else if (h > 0 && h < verNum - 1)
			return S1;
		else if (verNum - 1 == h)
			return S2;
		else
			return Error;
	}

	inline int NormalFormMatrix4HP::getLongestPath(){
		int i = 0;
		while (g_p[i][i + 1] > 0 && i < verNum - 1){
			i++;
		}
		return i+1;
	}

	inline void printPath(std::ofstream & of){
		of << "The current longest path is: ";
		printf("The current longest path is: ");
		int i = 0;
		double score = 1.0;
		while(g_p[i][i+1] > 0 && i < verNum - 1)
		{

			printf("%3i, ", v_heading[i]);
			of << v_heading[i] << ", ";
			score *= g_p[i][i + 1];
			i++;
		}
		printf("%3i.", v_heading[i]);
		of << v_heading[i] << ".";
		printf("\nIts score is %8.4f\n", score);
		of << "\nIts score is " << score << std::endl;
	}

	// before interchanging column h_i and h_j, 
	// and interchanging row h_i and row h_j;
	// four elements are involved along the normal path.
	// elements : (i -1, i), (i, i+1), (j -1, j), and (j, j+ 1);

	// after interchanging column h_i and h_j, 
	// and interchanging row h_i and row h_j;
	// four elements are involved along the normal path.
	// elements : (i -1, j), (j, i+1), (j -1, i), and (i, j +1);
	// **********
    // Return: 
	// **********
	//   the change of the non-zero elements before and 
	//   after the column/row interchanging.
	inline int delta_non_zero_elements(
		const int & h_i, // heading index i,
		const int & h_j, // heading index j,
		const int & e_im1_i, // matrix element at (i -1, i);
		const int & e_i_ip1, // matrix element at (i, i + 1);
		const int & e_jm1_j, // matrix element at (j -1, j);
		const int & e_j_jp1, // matrix element at (j, j +1);
		const int & e_im1_j, // matrix element at (i -1, j);
		const int & e_j_ip1, // matrix element at (j, i + 1);
		const int & e_jm1_i, // matrix element at (j -1, i);
		const int & e_i_jp1 // matrix element at (i, j +1);
		){
		
		int count1 = 0, count2 = 0;
		// for count1
		if (h_j - h_i == 1){ 
			// g_p[h_i][h_i + 1] <==> g_p[h_j - 1][h_j]
			count1 += e_im1_i == 0 ? 0 : 1;
			count1 += e_i_ip1 == 0 ? 0 : 1;
			count1 += e_j_jp1 == 0 ? 0 : 1;
		}
		else if (h_i - h_j == 1){ 
			// g_p[h_i - 1][h_i] <==> g_p[h_j][h_j + 1]
			count1 += e_im1_i == 0 ? 0 : 1;
			count1 += e_i_ip1 == 0 ? 0 : 1;
			count1 += e_jm1_j == 0 ? 0 : 1;
		}
		else{
			count1 += e_im1_i == 0 ? 0 : 1;
			count1 += e_i_ip1 == 0 ? 0 : 1;
			count1 += e_jm1_j == 0 ? 0 : 1;
			count1 += e_j_jp1 == 0 ? 0 : 1;
		} // end for count1

		// for count2
		count2 += e_im1_j == 0 ? 0 : 1;
		count2 += e_j_ip1 == 0 ? 0 : 1;
		count2 += e_jm1_i == 0 ? 0 : 1;
		count2 += e_i_jp1 == 0 ? 0 : 1;
		// end for count2
		return count2 - count1;
	}
};

#endif