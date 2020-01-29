/**
* @file: heurisHP.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "heurisHP.hpp"

/* This algorithm is proposed in this paper.
 * A HEURISTIC METHOD FOR THE DETERMINATION
 * OF A HAMILTONIAN CIRCUIT IN A GRAPH.
 * see http://journals.cambridge.org/download.php?file=%2FANZ%2FANZ28_03%2FS0334270000005439a.pdf&code=2fefdbb07e26bccf0e7fa9d615b280fd
 */

NormalFormMatrix4HP::NormalFormMatrix4HP(double ** g,
	const int & v) :
	verNum(v), g_p(g)
{
	v_out_neibor_distur =
		vector<vector<std::pair<int, int>>>(verNum);
	for (auto i = 0; i < verNum; ++i){
		for (auto j = 0; j < verNum; ++j){
			if (g_p[i][j] > 0)
				v_out_neibor_distur[i].push_back(
				make_pair(j, 0));
		}
	}
	// initialize the normal path of the graph matrix;
	for (auto i = 0; i < verNum - 1; ++i)
		v_normal_path.push_back(g_p[i][i + 1]);

	// initialize the index for row and column;
	for (auto i = 0; i < verNum; ++i)
		v_heading.push_back(i);
	for (auto i = 0; i < verNum; ++i)
		v_node.push_back(i);
}


bool NormalFormMatrix4HP::isNormalFormMat(int & c){
	for (int i = 0; i < verNum - 1; ++i){
		if (g_p[i][i + 1] == 0){
			c = i + 1;
			return false;
		}
	}
	c = -1;
	std::cout << "The matrix is already in normal form."
		<< " An HP exists.\n";
	return true;
}

void NormalFormMatrix4HP::swap_rows_clos
(const int & h_idx1, /*heading index 1*/
const int &h_idx2/*heading index 2*/)
{

	// interchange row idx1 and row idx2;
	// save row idx1 firstly.
	vector<double> v_temp(g_p[h_idx1], g_p[h_idx1] + verNum);
	// std::copy -- Copies the elements in the range[first, last) 
	// into the range beginning at result.
	// copy row idx2 to row idx1;
	for (auto i = 0; i < verNum; ++i){
		g_p[h_idx1][i] = g_p[h_idx2][i];
	}
	// copy row idx1 to row idx2;
	for (auto i = 0; i < verNum; ++i){
		g_p[h_idx2][i] = v_temp[i];
	}

	// interchange column idx1 and column idx2;
	// save column idx1 firstly.
	for (auto i = 0; i < verNum; ++i){
		v_temp[i] = g_p[i][h_idx1];
	}
	// copy column idx2 to column idx1;
	for (auto i = 0; i < verNum; ++i){
		g_p[i][h_idx1] = g_p[i][h_idx2];
	}
	// copy column idx1 to column idx2;
	for (auto i = 0; i < verNum; ++i){
		g_p[i][h_idx2] = v_temp[i];
	}

	// swap the order of the vertex
	// in v_heading;
	int temp_node = v_heading[h_idx1];
	v_heading[h_idx1] = v_heading[h_idx2];
	v_heading[h_idx2] = temp_node;
	// update the heading index
	// in v_node.
	v_node[v_heading[h_idx1]] = h_idx1;
	v_node[v_heading[h_idx2]] = h_idx2;
}

// print the data member matrix 
void NormalFormMatrix4HP::showMatrix(){
	for (int i = 0; i < verNum; ++i){
		printf("%6i, ", v_heading[i]);
	}
	std::cout << " (vertex idx)\n";
	for (int i = 0; i < verNum; ++i){
		for (int j = 0; j < verNum; ++j){
			printf("%6.4f, ", g_p[i][j]);
		}
		cout << std::endl;
	}
}

// print any matrix;
void NormalFormMatrix4HP::showMatrix(const double ** g_p_){
	for (int i = 0; i < verNum; ++i){
		for (int j = 0; j < verNum; ++j){
			printf("%6.4f, ", g_p_[i][j]);
		}
		cout << std::endl;
	}
}

// return how many non-zero nodes have been removed
// after interchanging row/column idx 
// and the candidate candiIdx.
// return value:
//  > 0 : the interchanging increases the normal path;
//  = 0 : the interchanging does not change the normal path;
//  > 0 : the interchanging decreases the normal path;
//   -1 : error occurs.
int NormalFormMatrix4HP::DisturbDegree(
	const int & h_i, /*the column/row heading h_i,
					 corresponding to the current zero element
					 encountered on the normal path.*/
					 const int & h_j) /*candidate heading h_j*/
{
	int count1 = 0, count2 = 0;
	// matrix elements involved in the normal path
	// when interchanging row/column h_i and h_j.
	// vector<int> v_nodes, v_new_nodes;
	// for h_i

	// the state if h_i / h_j:
	// 1) state 0: h_i / h_j = 0;
	// 2) state 1: h_i / h_j \in [1, verNum -2];
	// 3) state 2: h_i / h_j = verNum -1;
	// 4) state 3: error occurs.
	idx_state s_hi = getIdxState(h_i),
		s_hj = getIdxState(h_j);
	// 2-bit 3-based number (s_hi || s_hj).
	int sta = 3 * s_hi + s_hj;
	switch (sta){
	case 0: // s_hi = S0, s_hj = S0;
		// usually, h_i != h_j, so this case cannot often happen.
		return 0;
	case 1: // s_hi = S0, s_hj = S1;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			0,// matrix element at (i -1, i);
			g_p[h_i][h_i + 1],// matrix element at (i, i + 1);
			g_p[h_j - 1][h_j],// matrix element at (j -1, j);
			g_p[h_j][h_j + 1],// matrix element at (j, j +1);
			0,// matrix element at (i -1, j);
			g_p[h_j][h_i + 1],// matrix element at (j, i + 1);
			g_p[h_j - 1][h_i],// matrix element at (j -1, i);
			g_p[h_i][h_j + 1]// matrix element at (i, j +1);
			);
		break;
	case 2: // s_hi = S0, s_hj = S2;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			0,// matrix element at (i -1, i);
			g_p[h_i][h_i + 1],// matrix element at (i, i + 1);
			g_p[h_j - 1][h_j],// matrix element at (j -1, j);
			0,// matrix element at (j, j +1);
			0,// matrix element at (i -1, j);
			g_p[h_j][h_i + 1],// matrix element at (j, i + 1);
			g_p[h_j - 1][h_i],// matrix element at (j -1, i);
			0// matrix element at (i, j +1);
			);
		break;
	case 3: // s_hi = S1, s_hj = S0;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			g_p[h_i - 1][h_i],// matrix element at (i -1, i);
			g_p[h_i][h_i + 1],// matrix element at (i, i + 1);
			0,// matrix element at (j -1, j);
			g_p[h_j][h_j + 1],// matrix element at (j, j +1);
			g_p[h_i - 1][h_j],// matrix element at (i -1, j);
			g_p[h_j][h_i + 1],// matrix element at (j, i + 1);
			0,// matrix element at (j -1, i);
			g_p[h_i][h_j + 1]// matrix element at (i, j +1);
			);
		break;
	
	case 4: // s_hi = S1, s_hj = S1;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			g_p[h_i - 1][h_i],// matrix element at (i -1, i);
			g_p[h_i][h_i + 1],// matrix element at (i, i + 1);
			g_p[h_j - 1][h_j],// matrix element at (j -1, j);
			g_p[h_j][h_j + 1],// matrix element at (j, j +1);
			g_p[h_i - 1][h_j],// matrix element at (i -1, j);
			g_p[h_j][h_i + 1],// matrix element at (j, i + 1);
			g_p[h_j - 1][h_i],// matrix element at (j -1, i);
			g_p[h_i][h_j + 1]// matrix element at (i, j +1);
			);
		break;
	
	case 5: // s_hi = S1, s_hj = S2;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			g_p[h_i - 1][h_i],// matrix element at (i -1, i);
			g_p[h_i][h_i + 1],// matrix element at (i, i + 1);
			g_p[h_j - 1][h_j],// matrix element at (j -1, j);
			0,// matrix element at (j, j +1);
			g_p[h_i - 1][h_j],// matrix element at (i -1, j);
			g_p[h_j][h_i + 1],// matrix element at (j, i + 1);
			g_p[h_j - 1][h_i],// matrix element at (j -1, i);
			0// matrix element at (i, j +1);
			);
		break;
	case 6: // s_hi = S2, s_hj = S0;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			g_p[h_i - 1][h_i],// matrix element at (i -1, i);
			0,// matrix element at (i, i + 1);
			0,// matrix element at (j -1, j);
			g_p[h_j][h_j + 1],// matrix element at (j, j +1);
			g_p[h_i - 1][h_j],// matrix element at (i -1, j);
			0,// matrix element at (j, i + 1);
			0,// matrix element at (j -1, i);
			g_p[h_i][h_j + 1]// matrix element at (i, j +1);
			);
		break;
	case 7: // s_hi = S2, s_hj = S1;
		return delta_non_zero_elements(
			h_i, // heading index i,
			h_j, // heading index j,
			g_p[h_i - 1][h_i],// matrix element at (i -1, i);
			0,// matrix element at (i, i + 1);
			g_p[h_j - 1][h_j],// matrix element at (j -1, j);
			g_p[h_j][h_j + 1],// matrix element at (j, j +1);
			g_p[h_i - 1][h_j],// matrix element at (i -1, j);
			0,// matrix element at (j, i + 1);
			g_p[h_j - 1][h_i],// matrix element at (j -1, i);
			g_p[h_i][h_j + 1]// matrix element at (i, j +1);
			);
		break;
	case 8: // s_hi = S2, s_hj = S2;
		return 0;
	default:{
		cout << "Error in function DisturbDegree().\n";
		return -1;
	}
	}
}

/*the column/row heading h_i,
corresponding to the current zero element
encountered on the normal path.
Return: the node index with the least disturbed degree.
*/
int NormalFormMatrix4HP::least_disturbed_node
(const int & h_i)
{
	// Notation:
	//   v_heading[h_i] = vertex_idx;
	// row-heading index = column-heading index - 1;
	// therefore we use h_i - 1 here.
	int nodeIdx = v_heading[h_i - 1];
	// the node with least disturbed degree.
	int least_disturb_nodeIdx = verNum + Large_Value,
		max_dist_degree = (-1) * Large_Value,
		temp_dist_degree = 0;
	for (auto temp_node_dist :
		v_out_neibor_distur[nodeIdx]){

		// /* int temp_neibor = temp_node_dist.first;
		//  * int temp_dist   = temp_node_dist.second;
		//  * int temp_h_i    = v_node[temp_neibor];
		//  */
		temp_dist_degree =
			DisturbDegree(h_i, v_node[temp_node_dist.first]);
		// update disturb degree;
		temp_node_dist.second = temp_dist_degree;
		if (max_dist_degree <= temp_dist_degree){
			max_dist_degree = temp_dist_degree;
			least_disturb_nodeIdx = temp_node_dist.first;
		}
	}
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	cout << " For Heading Idx " << h_i
		<< ", and the corresponding vertex "
		// row-heading index = column-heading index - 1;
		// therefore we use h_i - 1 here.
		<< v_heading[h_i - 1]
		<< ", the least disturb vertex is "
		<< least_disturb_nodeIdx << ", with disturb degree = "
		<< max_dist_degree << endl;
#endif
	return least_disturb_nodeIdx;
}