/**
* @file: greedyDFS.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_Greedy_DFS_H_
#define __HEADER__CROWD_SOURCING_Greedy_DFS_H_
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <string>
#include <stack>          // std::stack
#include "path.hpp"
using namespace std;
//typedef std::pair<double, Path> ScorePathPAIR;
typedef std::pair<double, vector<int>> ScorePathPAIR;


class Greedy_DFS_Path{
private:
	// In DFS, each vertex has three possible 
	// colors representing its state:
	//   * white : vertex is unvisited;
	//   * gray: vertex is in progress;
	//   * black: DFS has finished processing the vertex.
	const int V; // vertex number
	const bool isDisplay;
	const bool isSaveData;
	double ** p_Gp;
	// Pointer to an array containing adjacency lists
	std::list<std::pair<int, double>> *adj;
	enum VertexColors : int { BLACK = 1, WHITE = 2, GRAY = 3 };
	vector<VertexColors> v_Colors;
	vector<ScorePathPAIR> v_HP;
	std::stack<int> s_tempPath; // temporal stack path
	std::stack<int> s_longestPath; // the possible longest path if no HP;
	vector<int> v_tempPath; // temporal vector path
	int temp_max_path_leng = 0;
	bool isHP = false;
public:
	Greedy_DFS_Path(double * * g_,
		const int & vNum, // vertex number;
		const bool & isDisplay_,
		const bool & isSaveData_);
	Greedy_DFS_Path::~Greedy_DFS_Path();
	void PathDFS_version2(const int & s, std::ofstream  & of);
	void PathDFS_version1(const int & s, std::ofstream  & of);
	void setlabel(const int & u, const VertexColors & color);
	void sortAdjListWeights();

	// for all the vertices v in V
	// search HPs beginning with each vertex;
	void GreedyPathDFS(std::ofstream  & of);

	// for the vertex v labeled by the input vertex index;
    // this function is used to run in parallel.
	// search an HP beginning with some specific vertex;
	void GreedyPathDFS(const int & v_idx, std::ofstream  & of);
	void printHPs(std::ofstream  & of);
	template < typename T > void printStack(const std::stack<T>& stk, std::ofstream  & of);
	template < typename T > void savePrintStack(const std::stack<T>& stk,
		std::ofstream  & of);
	//std::ofstream  of;
	double getScore();
}; /*end of class*/

	// descending sorting
	bool desSort(const pair<int, double> & p1,
		const pair<int, double> & p2);
	

#endif