/**
* @file: path.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_PROJECT_PATH_H_
#define __HEADER__CROWD_SOURCING_PROJECT_PATH_H_
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
#include "myUtility.hpp" /*factorial()*/

// #define K 20 // top-k ranking, we k = 1 currently
// #define VertexNUM 20
// path length means how many edges it contains;
// #define PathLength VertexNUM -1 

typedef unsigned int UINT;
typedef unsigned int uint;
typedef unsigned long WEIGHT_ID;
typedef std::pair<unsigned long, double> id_weights_pair_type;
typedef std::vector<std::pair<unsigned long, double>>::iterator 
v_weights_iterator_type;
typedef unsigned long ulong;
typedef unsigned long long ullong;
using namespace std;

struct Path { /*struct Path definition*/
	 unsigned int vertexNum;
	 unsigned int pathLen;
	// unsigned int pathID;
	std::vector <int> v_nodes;
	double score = 1.0;

	Path(const UINT & vNum) : 
		vertexNum(vNum), pathLen(vNum - 1){
		score = 1.0;
	}

	/* No assignment (=) operator is defined for the class. 
	 * If you define any assignment operator that takes the 
	 * class as a parameter, the compiler cannot generate a 
	 * default assignment operator for that class.
     * Assignment operators are not inherited by derived classes. 
	 * You must explicitly define an assignment operator 
	 * for each class.
	 */
	Path & Path:: operator= (const Path & pSource);

	// Copy constructor
	Path(const Path & p): 
	vertexNum(p.vertexNum),
	pathLen (p.pathLen),
	score (p.score){
		v_nodes = std::vector<int>(p.v_nodes);
	}

	Path(const UINT & pLen, const UINT & vNum):
		vertexNum(vNum),
	    pathLen (pLen) {}

	~Path(){
		// vector::swap can be used to 
		// effectively free allocated memory.
		vector<int >().swap(v_nodes);
	};

	// 2D array of weights;
	inline  void  setScore(double & s){
		score = s;
	};

	inline double getScore(){ return score; }
	void vertex2Path(int * vertices);
	void printPath(ullong & pathID);
	void Path::printPath(std::ofstream  & of);
}; /*end of struct Path definition*/


typedef std::pair<double, Path> score_path_pair_type;

bool isLarger(const id_weights_pair_type & arg1,
	const id_weights_pair_type & arg2);

template <typename T1, typename T2>
bool isLess(const std::pair<T1, T2> & i,
	const std::pair<T1, T2> & j){
	return (i.first < j.first);
}

//bool isLess(const std::pair<double, Path> &, const std::pair<double, Path> &);
#endif