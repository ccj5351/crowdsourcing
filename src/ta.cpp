/**
* @file: ta.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#include "ta.hpp"
using namespace std;

bool TA::getScoreWithThre(std::vector<int>  & tempVertexSet,
	const double & thre, double & temScore){
	temScore = 1.0;
	auto tsize = tempVertexSet.size() - 1;
	int i = 1;
	for (vector<int>::iterator it = tempVertexSet.begin();
		it != tempVertexSet.end() - 1; ++it, i++){
		temScore *= v_weights[(*it) * vertexNum + (*(it + 1))];
		if (temScore <= thre){
			temScore = thre;
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			if (IS_DISPLAY_STH){
				std::cout << " >- Early stop at edge " << i << "/" << tsize
					<< ", for the currently scanning"
					<< " path, already smaller than the minimum"
					<< " score among those top " << top_k << " paths.\n";
			}
#endif
			return false;
		}
	}
	return true;
}

bool TA::isPath(
	// complete vertex set, if some vertex has been
	// removed, it should be repaired and added.
	V_M_TYPE & RepairedVertexSet)
{
	for (V_M_TYPE::iterator it = RepairedVertexSet.begin();
		it != RepairedVertexSet.end() - 1; ++it){
		if (v_weights[(*it) * vertexNum + (*(it + 1))] <= 0){
            #ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			if (IS_DISPLAY_STH)
				cout << "This vertex set permutation cannot make a path,"
				<< " since no edge exists from " << (*it) << " to "
				<< *(it + 1) << endl;
            #endif
			return false;
		}
	}
	return true;
}


double TA::getScoreWOThre(std::vector<int> & tempVertexSet)
{
	double temScore = 1.0;
	for (vector<int>::iterator it = tempVertexSet.begin();
		it != tempVertexSet.end() - 1; ++it){
		temScore *= v_weights[(*it) * vertexNum + (*(it + 1))];
	}
	return temScore;
}

// sort the weights in the descending order;
void TA::getSortedWeights(vector<double> & v_wgts){
	int rIdx = 0, cIdx = 0;
	for (int w_id = 0, s = v_wgts.size();
		w_id < s; ++w_id){
		rIdx = w_id / vertexNum; // row index
		cIdx = w_id % vertexNum; // column index
		v_sorted_weights.push_back(
			make_pair(w_id, v_wgts[w_id]));
	}


	// sort the list of weights;
	// sort by value using std::sort,
	// complexity is O(n*log(n))
	bool(*fn_pt)(const id_weights_pair_type &,
		const id_weights_pair_type &) = isLarger;
	std::sort(v_sorted_weights.begin(),
		v_sorted_weights.end(), fn_pt);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	if (IS_DEBUG_MODE){
		cout << "v_sorted_weights = \n";
		for (vector<id_weights_pair_type>::iterator
			it = v_sorted_weights.begin();
			it != v_sorted_weights.end(); ++it){
			std::cout << it->first << " => "
				<< it->second << '\n';
		}
	}
#endif
}

void TA::initiaWeights(vector<double> & v_weights){
	this->v_weights = vector<double>(v_weights);
}

void TA::initiaWeights(double ** mat){
	if (v_weights.empty())
		v_weights = vector<double>(vertexNum*vertexNum);
	for (int i = 0; i < vertexNum; ++i)
		for (int j = 0; j < vertexNum; ++j)
			this->v_weights[i*vertexNum + j] = mat[i][j];
}

// do allocating and initialization 
// for the top-k paths; 
void TA::initiTopKPaths(){
	Path tp(pathLen, vertexNum);
	for (uint i = 0; i < top_k; i++){
		v_top_k_paths.push_back(std::pair<double, Path>(0, tp));
	}
}

void TA::bruteForceParallel(
	const ullong & left, // staring point(including);
	const ullong & right, // ending point(including);
	std::ofstream  & of){
	if (!of){
		cout << "failed open file " << endl;
	}
	string l_combi, r_combi;
	double min_topK_score = Large_Value;
	Path temPath(pathLen, vertexNum);
	ullong idx = 0, TN = (right - left),
		percent1 = TN / vertexNum,
		temp_permut = left;
	Fact_Int_Mapping f2i(vertexNum, true);
	std::cout << "The staring point:\n";
	of << "The staring point:\n";
	// starting from the left point;
	vector<int> v_temp_combi = f2i.int2fact2(left, l_combi, of);
	std::cout << "\nThe ending point:\n";
	of << "\nThe ending point:\n";
	vector<int> v_right_combi = f2i.int2fact2(right, r_combi, of);
	std::cout << "\nDoing path searching...\n";
	of << "\nDoing path searching...\n";
	// sort in ascending order
	bool(*fn_pt)(const std::pair<double, Path> &,
		const std::pair<double, Path> &) = isLess < double, Path > ;
	do
	{ // permutation of vertexSet 

		if ((idx + 1) % percent1 == 0){
			f2i.isDisplay = true;
			std::cout << "****************\nFinished "
				<< (idx + 1) / percent1
				<< "/" << vertexNum << "\n";
			of << "****************\nFinished "
				<< (idx + 1) / percent1
				<< "/" << vertexNum << "\n";
			ullong re = f2i.fact2int(&v_temp_combi[0], of);
			f2i.int2fact2(re, l_combi, of);
		}
		idx++;
		temp_permut++;
		// if this permutation of vertex set makes a path,
		// that is, there are not any zero-weight edges 
		// in this permutation.
		if (isPath(v_temp_combi)){ // be a path;
			temPath.vertex2Path(&v_temp_combi[0]);
			temPath.score = getScoreWOThre(v_temp_combi);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_			
			if (IS_DISPLAY_STH){
				temPath.printPath(idx);
			}
#endif

			// save the first top-k paths;
			if (v_top_k_paths.size() < top_k){
				v_top_k_paths.push_back(make_pair(temPath.score, temPath));
				min_topK_score = std::min(min_topK_score, temPath.score);
			}
			else {
				// update the paths
				if (min_topK_score < temPath.score){
					// sort the top-k paths
					// Sorts the elements in the range[first, last) 
					// into ascending order.
					std::sort(v_top_k_paths.begin(), v_top_k_paths.end(), fn_pt);
					v_top_k_paths[0] = make_pair(temPath.score, temPath);
					// update the min value of top-k paths;
					min_topK_score = std::min(temPath.score,
						v_top_k_paths[1].first);
				}

			}
		}
		else{ // no path for this permutation of vertex set.

			// The "continue" statement causes the program to skip 
			// the rest of the loop in the current iteration, 
			// as if the end of the statement block had been reached, 
			// causing it to jump to the start of the following iteration.
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			if (IS_DISPLAY_STH)
				std::cout << "This permutation of vertex set cannot make a path"
				<< ", 'continue' jumps to next round of While-Loop.\n";
#endif
			continue;
		}

	} // end of while
	while (std::next_permutation(v_temp_combi.begin(), v_temp_combi.end(),
	// pay attention to this function
	// starting from ascending order, 
	// for example, begin with "0, 1, 2, 3 ,4",
	// ending up with descending order, i.e., "4, 3, 2, 1, 0".
	std::less<int>()) && (temp_permut <= right));
	std::sort(v_top_k_paths.begin(), v_top_k_paths.end(), fn_pt);
	if (v_top_k_paths.empty())
		of << "No HP exists.\n";
	else{
		of << "HP is :\n";
		for (int i = 0; i < v_top_k_paths.size(); ++i)
			v_top_k_paths[i].second.printPath(of);
	}

}
int * TA::getmaxPath(){
	int * pt = new int[vertexNum];
	for (int i = 0; i < vertexNum; ++i)
		pt[i] = MaxPath->v_nodes[i];
	return pt;
}

// save the top-1 path;
bool TA::bruteForce(std::vector<int> & vertexSet,
	std::ofstream  & of){
	if (!of){
		cout << "failed open file " << endl;
	}
	double max_score = -1;
	ullong path_counter = 0;
	Path temPath(pathLen, vertexNum);
	ullong TN = factorial(vertexNum);
	isHP = false;
	do
	{ // permutation of vertexSet 
		// if this permutation of vertex set makes a path,
		// that is, there are not any zero-weight edges 
		// in this permutation.
		path_counter++;

		if (path_counter == TN * 0.1)
			std::cout << "Finish 10%\n";
		if (path_counter == TN * 0.2)
			std::cout << "Finish 20%\n";
		if (path_counter == TN * 0.3)
			std::cout << "Finish 30%\n";
		if (path_counter == TN * 0.4)
			std::cout << "Finish 40%\n";
		if (path_counter == TN * 0.5)
			std::cout << "Finish 50%\n";
		if (path_counter == TN * 0.6)
			std::cout << "Finish 60%\n";
		if (path_counter == TN * 0.7)
			std::cout << "Finish 70%\n";
		if (path_counter == TN * 0.8)
			std::cout << "Finish 80%\n";
		if (path_counter == TN * 0.9)
			std::cout << "Finish 90%\n";


		if (isPath(vertexSet)){ // be a path;
			isHP = true;
			temPath.vertex2Path(&vertexSet[0]);
			temPath.score = getScoreWOThre(vertexSet);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			if (IS_DISPLAY_STH){
				temPath.printPath(path_counter);
			}
#endif
			if (max_score < temPath.score){
				max_score = temPath.score;
				*MaxPath = temPath;
			}
		}
		else{ // no path for this permutation of vertex set.

			// The "continue" statement causes the program to skip 
			// the rest of the loop in the current iteration, 
			// as if the end of the statement block had been reached, 
			// causing it to jump to the start of the following iteration.
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			if (IS_DISPLAY_STH)
				std::cout << "This permutation of vertex set cannot make a path"
				<< ", 'continue' jumps to next round of While-Loop.\n";
#endif
			continue;
		}
		
	} // end of while
	while (std::next_permutation(vertexSet.begin(), vertexSet.end(),
	// pay attention to this function
	// starting from ascending order, 
	// for example, begin with "0, 1, 2, 3 ,4",
	// ending up with descending order, i.e., "4, 3, 2, 1, 0".
	std::less<int>()));
	if (isHP){
	std::cout << "Max score = " << max_score << endl;
	MaxPath->printPath(of);
	}
	else{
		cout << "No HP!\n";
		of << "No HP!\n";
	}
	return isHP;
}

ullong TA::getTopKPaths(
	std::vector<int> & vertexSet) // vertex set;
{

	initiTopKPaths();
	// scanning the sorted weights of edges in descending order;
	double thre = .0; // threshold value;
	uint weight_row_id = 0, // parent node id;
		weight_col_id = 0; // child node id;
	// thus e.g., (weight_row_id, weight_col_id) = (2, 3) 
	// means the directed edge 
	// form vertex 2 to vertex 3 (0-based indexing).
	Path temPath(pathLen, vertexNum);
	ullong path_counter = 0;
	ullong top_k_idx = 0;
	// sort in ascending order
	bool(*fn_pt)(const std::pair<double, Path> &,
		const std::pair<double, Path> &) = isLess<double,Path>;

	// auxiliary vertex set
	int edge_idx = 1,
		edges_size = v_sorted_weights.size();
	min_of_top_K_Scores = -1 * Large_Value;
	for (v_weights_iterator_type it = v_sorted_weights.begin();
		it != v_sorted_weights.end(); ++it, ++edge_idx)
	{ /*current scanning */
		w_id = it->first;

		// for example, currently scanning the weight of the 
		// edge (vertex 2, vertex 3)
		// weight = w(2, 3), thus we should enumerate all the 
		// paths involving the
		// edge of (v2, v3), which can be obtained by permuting 
		// all the vertices, but v2 and v3 should be bound together
		// during the permutation.
		weight_row_id = w_id / vertexNum; // e.g., = 2
		weight_col_id = w_id % vertexNum; // e.g., = 3

		// How to bind vertex 2 and vertex 3?
		// We can just remove vertex 3 from the vertex set, 
		// but keep vertex 2, then do the permutation,
		// in order to get all the possible results.
		// Then for each path, please add vertex 3 closely 
		// after vertex 2. In other words, here vertex 2 actually
		// stands for the combination of vertex 2 and vertex 3.

		// "bind" parent vertex and child vertex by removing 
		// the child vertex;
		/*  erase the 6th element, that it erase the element
		 *  with the index 5, zero-based indices.
		 *  myvector.erase (myvector.begin()+5);
		 */
		/*
		 * for example, weight_col_id = 3;
		 * this operation will erase the 4-th element in the vector.
		 * and the 4-th element is exactly weight_col_id (= 3) here.
		 * Due to the predefined relationship between the index and the element,
		 * vertexSet[0] = 0 (i.e., vertex 0),
		 * vertexSet[1] = 1 (i.e., vertex 1),
		 * vertexSet[2] = 2 (i.e., vertex 3),
		 * vertexSet[3] = 3 (i.e., vertex 4), etc.
		 */
		// erase operation's complexity is linear 
		// on the number of elements erased (destructions) plus
		// the number of elements after the last element deleted (moving).
		std::vector<int> vertexSetCopy(vertexSet);
		vertexSet.erase(vertexSet.begin() + weight_col_id);
		thre = pow(it->second, vertexNum - 1);
		cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
		std::cout << "Currently " << edge_idx << "/" << edges_size
			<< " scanning edge (" << weight_row_id
			<< ", " << weight_col_id << "), weight = "
			<< it->second << ", Threshold = " << thre << endl;

		// after erase one element, do sorting for 
		// the following permutation.
		// Sorts the elements in the range[first, last) 
		// into ascending order.
		std::sort(vertexSet.begin(), vertexSet.end());
		// do-while will first print array values, then do next permutation;
		// but while-do will firstly do next permutation, 
		// then print array values.

		do{ // permutation of vertexSet 

			V_M_TYPE  repairedVertexSet
				= repairVertexSet(vertexSet, weight_row_id,
				weight_col_id);

			// TA need stop or not
			if (thre < min_of_top_K_Scores){
				cout << "TA success!\n";
				for (uint i = 0; i < top_k; ++i){
					cout << v_top_k_paths[i].second.score << ", ";
				}
				cout << endl;
				return path_counter;
			}

			// if this permutation of vertex set makes a path,
			// that is, there are not any zero-weight edges 
			// in this permutation.
			if (isPath(repairedVertexSet)){ // be a path;
				/* save the top-k paths
				* firstly, to save the first generated k paths;
				*/
				if (top_k_idx < top_k){
					double temScore = getScoreWOThre(repairedVertexSet);
					temPath.setScore(temScore);
					temPath.vertex2Path(&repairedVertexSet[0]);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
					if (IS_DISPLAY_STH){
						temPath.printPath(path_counter);
					}
#endif
					// save the first k paths;
					v_top_k_paths[top_k_idx].first = temScore;
					// error C2582: 'operator =' function is unavailable in ...
					// v_top_k_paths[Idx_min_of_K_Scores].second = Path(temPath);
					v_top_k_paths[top_k_idx].second = temPath;
				}
				else { // for the newly generated path 
					   // after the first k paths.
					double temScore = 1.0;
					// then, compare the newly generated path, with the 
					// minimum score of top-k paths.

					//std::pair<int, double> temp = min_of_K_score_path_pair();
					//int idx_min_k_scores = temp.first;
					//min_of_K_Scores = temp.second;

					if (getScoreWithThre(repairedVertexSet,
						min_of_top_K_Scores, temScore)){
						// if TRUE, means threshold actually does not work,
						// i.e., "threshold (= min_of_K_Scores here) < temScore;
						// Thus, save the value temScore.
						temPath.setScore(temScore);
						temPath.vertex2Path(&repairedVertexSet[0]);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
						if (IS_DISPLAY_STH){
							temPath.printPath(path_counter);
						}
#endif
						// sort the top-k paths
						// Sorts the elements in the range[first, last) 
						// into ascending order.
						std::sort(v_top_k_paths.begin(), v_top_k_paths.end(), 
							fn_pt);
						// Thus, the first element has the minimum score.
						min_of_top_K_Scores = v_top_k_paths[0].first;

						// save the currently scanned path
						// since its score is larger than 
						// the minimum one among
						// the top-k paths.
						v_top_k_paths[0].first = temScore;

						// error C2582: 'operator =' function 
						// is unavailable in ...
						// v_top_k_paths[Idx_min_of_K_Scores].second 
						// = Path(temPath);
						v_top_k_paths[0].second = temPath;

						// update the min value of top-k paths;
						// since maybe the element at index [1] 
						// is larger than temScore;
						min_of_top_K_Scores = std::min(temScore,
							v_top_k_paths[1].first);
					}
				} /*end of updating top-k paths*/
				top_k_idx++;
				path_counter++;
		} /* end of being a path;*/
			
			else{ // no path for this permutation of vertex set.

				// The "continue" statement causes the program to skip 
				// the rest of the loop in the current iteration, 
				// as if the end of the statement block had been reached, 
				// causing it to jump to the start of the following iteration.
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
				if (IS_DISPLAY_STH)
					cout << "This permutation of vertex set cannot make a path"
					<< ", 'continue' jumps to next round of While-Loop.\n";
#endif
				continue;
			}
		} while (std::next_permutation(vertexSet.begin(), vertexSet.end(),
			// pay attention to this function
			// starting from ascending order, 
			// for example, begin with "0, 1, 2, 3 ,4",
			// ending up with descending order, i.e., "4, 3, 2, 1, 0".
			std::less<int>()));

		// for next round scanning, repair the vertex 
		// set to the original set.
		vertexSet = std::vector<int>(vertexSetCopy);

		// TA need stop or not
		if (thre < min_of_top_K_Scores){
			cout << "TA success!\n";
			for (uint i = 0; i < top_k; ++i){
				cout << v_top_k_paths[i].second.score << ", ";
			}
			cout << endl;
			return path_counter;
		}


	} /*end of for, i.e., current scanning*/

	return path_counter;
};

bool TA::savePath(std::ofstream  & of){
	if (!of){
		cout << "failed open file " << endl;
		return false;
	}
	else{
		if (v_top_k_paths.empty()){
			of << "No HP exists.\n";
		}
		else {
			of << "Top " << top_k << " paths include:\n";
			for (int i = 0; i < top_k; ++i){
				for (uint j = 0; j < vertexNum - 1; j++)
					of << v_top_k_paths[i].second.v_nodes[j] << " --> ";
				of << v_top_k_paths[i].second.v_nodes[vertexNum - 1]
					<< ", score = " << v_top_k_paths[i].first << "\n";
			}
		}
		return true;
	}
}

/*
// example showing that how to use this function.
int main(){
double P[] =
{ 1.0, 0.8, 0.8, 0.7273,
0.2, 1.0, 0.5, 0.4,
0.2, 0.5, 1.0, 0.4,
0.2727, 0.6, 0.6, 1.0
};
bool isDisplay = false;
bool is_Debug = false;
std::vector<unsigned int> vertexSet = { 0, 1, 2, 3 };
int VertexNUM = vertexSet.size();
std::vector<double> v_weights(P, P + VertexNUM * VertexNUM);
int K = 50;
// initialize static variables in class TA;
TA ta(VertexNUM, K, isDisplay, is_Debug);
ta.initiTopKPaths(v_weights);
ta.getTopKPaths(vertexSet, v_weights);
return 0;
}
*/