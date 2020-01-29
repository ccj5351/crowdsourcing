/**
* @file: ta.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_PROJECT_TA_H_
#define __HEADER__CROWD_SOURCING_PROJECT_TA_H_
#include <iostream>
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
#include <fstream>      // std::ofstream
#include "path.hpp"
#include "factorial2int.hpp"

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_

//#define K 20 // top-k ranking, we k = 1 currently
//#define VertexNUM 20
//#define PathLength VertexNUM -1 // path length means how many edges it contains;
//#define IS_DISPLAY_STH false // for display something
#define Large_Value 200 /*used to initialize TA::min_of_K_Scores*/
//#define IS_DEBUG_MODE true // for debugging
using namespace std;
typedef std::pair<double, Path> score_path_pair_type;
typedef std::vector <int> V_M_TYPE; // vertex label type.
class TA {// Threshold algorithm
private:
	const uint  vertexNum;
	const uint  top_k;
	const uint  edgeNum;
	const uint  pathLen;
	const ullong MaxPathNum;
	// for display something
	const bool IS_DISPLAY_STH;
	// for debugging
	const bool IS_DEBUG_MODE;
	//uint  Idx_min_of_K_Scores = 0;
	// whatever as long as it's larger than 1;
	double min_of_top_K_Scores = .0; 
	bool isHP;

	// const unsigned long duplicate_weights_num;
	WEIGHT_ID w_id = 0;

	// std::map V.S. std::unordered_map
	/*
	1. std::map containers
	Internally, the elements in a map are always sorted by 
	its key following a specific strict weak ordering criterion
	indicated by its internal comparison object (of type Compare).
	map containers are generally slower than unordered_map 
	containers to access individual elements by their key, 
	but they allow the direct iteration on subsets based on 
	their order.
	2. std::unordered_map containers
	unordered_map containers are faster than map containers to access
	individual elements by their key,although they are generally less
	efficient for range iteration through a subset of their elements.
	*/

	//std::unordered_map<unsigned long, double> m_sorted_weights;
	vector<std::pair<unsigned long, double>>  v_sorted_weights;
	// original weights matrix, i.e., unsorted;
	vector<double>  v_weights;
	std::vector<std::pair<double, Path>>  v_top_k_paths;
	Path * MaxPath = NULL;

public:
	inline ulong getPathNum(){ return MaxPathNum;}
	// first : path index within [0, top_k -1];
	// second : path score;
	inline  pair<int,double> min_of_K_score_path_pair(){
		int min_k_Scores = Large_Value;
		int Idx_min_of_K_Scores_ = 0;
		for (uint i = 0; i < top_k; ++i){
			if (min_k_Scores >= v_top_k_paths[i].first){
				min_k_Scores = v_top_k_paths[i].first;
				Idx_min_of_K_Scores_ = i;
			}
		}
		return make_pair(Idx_min_of_K_Scores_, min_k_Scores);
	}

	// constructor
	TA ( const uint & vNum, 
		const uint & e_, 
		const uint & k, 
		const bool & is_display,
	    const bool & id_debug) 
		: vertexNum(vNum), 
		edgeNum(e_),
		top_k(k),
		pathLen(vNum - 1),
		MaxPathNum(factorial(vNum)),
		IS_DISPLAY_STH(is_display),
		IS_DEBUG_MODE(id_debug){
		MaxPath = new Path(pathLen, vertexNum);
	}

	// combine column-index vertex and row-index vertex together.
	// specifically, append column-index vertex to row-index vertex,
	// making the edge (row-index, col-index).
	inline std::vector<int> repairVertexSet
		( std::vector<int> & vertexSet,
		const uint & rIdx,  // row-index vertex
		const uint & cIdx // column-index vertex.
		){

		std::vector<int>  tempVertexSet;
		for (vector<int>::iterator it = vertexSet.begin();
			it != vertexSet.end(); ++it){
			// for tempVertexSet
			tempVertexSet.push_back(*it);
			if ((*it) == rIdx){
				tempVertexSet.push_back(cIdx);
			}
		}
		return tempVertexSet;
	}

	bool isPath(
		// complete vertex set, if some vertex has been
		// removed, it should be repaired and added.
		V_M_TYPE & RepairedVertexSet
		);

	bool getScoreWithThre(std::vector<int>  & tempVertexSet,
		const double & thre, double & temScore);

	double getScoreWOThre(std::vector<int>  & tempVertexSet);

	void getSortedWeights(vector<double> & v_wgts);
	void initiTopKPaths();
	void initiaWeights(vector<double> & v_weights);
	void initiaWeights(double ** mat);
	bool bruteForce(std::vector<int> & vertexSet,
		std::ofstream  & of);
	int * TA::getmaxPath();
	void bruteForceParallel(
		const ullong & left, // staring point(including);
		const ullong & right, // ending point(including);
		std::ofstream  & of);

	ullong getTopKPaths(std::vector<int> & vertexSet);

	bool savePath(std::ofstream  & of);
};
#endif