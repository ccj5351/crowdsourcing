/**
* @file: myGraph.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_GRAPH_H_
#define __HEADER__CROWD_SOURCING_GRAPH_H_
// Below is a implementation of Graph Data Structure 
// in C++ as Adjacency List.
// http://stackoverflow.com/questions/5493474/graph-implementation-c
// https://www.cs.princeton.edu/~rs/Algs3.cxx5/code.txt
/*
* This file contains the code from "Algorithms in C++, Third Edition,
* Part 5," by Robert Sedgewick, and is covered under the copyright
* and warranty notices in that book. Permission is granted for this
* code to be used for educational purposes in association with the text,
* and for other uses not covered by copyright laws, provided that
* the following notice is included with the code:
* "This code is from "Algorithms in C++, Third Edition,"
* by Robert Sedgewick, Addison-Wesley, 2002."
* Commercial uses of this code require the explicit written
* permission of the publisher. Send your request for permission,
* stating clearly what code you would like to use, and in what
* specific way, to: aw.cse@aw.com
* */

#include <cstdio>      /* printf, scanf, puts, NULL */
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <iostream>
#include <climits>
#include <cstring> /* memset */
#include <string>
#include <queue>
#include <set>
#include <vector>
#include <fstream>      // std::fstream
#include <random>  /*std::uniform_int_distribution*/
#include <map> /*ordered map*/
#include <chrono> //std::chrono::system_clock

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_
using namespace std;
//using namespace Eigen;

// color 
// In DFS, each vertex has three possible 
// colors representing its state:
//   * white : vertex is unvisited;
//   * gray: vertex is in progress;
//   * black: DFS has finished processing the vertex.
enum VertexColors : int { BLACK = 1, WHITE = 2, GRAY = 3 };

struct vertex{
	typedef pair<int, vertex*> ve;
	vector<ve> adj; //cost of edge, destination vertex
	std::string name;
	// constructor
	vertex(string s){
		name = s;
	}
};

// 
struct Edge{
	int u, v;
	Edge(int u = -1, int v = -1) : u(u), v(v) { }
};


class myGraph{
private:
	typedef map<string, vertex *> vmap;
	vmap work;
	vector<int> v_degree;
	const int vertexNum;
	const bool isDirected = false;
	//int edgeNum;
	int degree_min;
	int degree_max;
	// In C++11 it is possible:
	// auto array = new double[M][N];
	// This way, the memory is not initialized. 
	// To initialize it do this instead:
	bool **  matrix_Gt = NULL;
	//vector<Edge> v_undirectedEdges;
	std::map<int, Edge> m_undirectedEdges;
	// Create a visited array and mark all vertices as not visited
	vector<VertexColors> v_Colors;
	vector<int> v_parents;
	vector<int> tempPath;
	const int path_Vertex_Length = 5;
	const int cycle_Vertex_Length = 4;
	int path_Edge_length = path_Vertex_Length - 1;
	std::vector<int> temp_vertexSet;

public:
	myGraph(int vNum, bool isdirect, int d_min, int d_max,
		int path_Vertex_Length_, int cycle_Vertex_Length_)
		: vertexNum(vNum), isDirected(isdirect),
		degree_max(d_max), degree_min(d_min),
		path_Vertex_Length(path_Vertex_Length_),
		cycle_Vertex_Length(cycle_Vertex_Length_){
		v_degree = vector<int>(vNum, 0);
		v_parents = vector<int>(vNum, -1);
		v_Colors = vector<VertexColors>(vNum, WHITE);
		matrix_Gt = new bool *[vertexNum];
		for (int i = 0; i < vertexNum; ++i){
			matrix_Gt[i] = new bool[vertexNum];
		}

		// initialize
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				*(matrix_Gt[i] + j) = false;
			}
		}
	}

	~myGraph(){
		for (int i = 0; i < vertexNum; ++i){
			delete[] matrix_Gt[i];
		}

		delete[] matrix_Gt;
		vector<int>().swap(v_degree);
		vector<int>().swap(v_parents);
		vector<VertexColors>().swap(v_Colors);
		vector<int>().swap(tempPath);
	}
	void setSink(const int & v);
	void setSource(const int & v);
	void setIsolateVertex(const int & v);
	int V() const { return vertexNum; }
	int E() const { return m_undirectedEdges.size(); }
	bool IsDirect() const { return isDirected; }
	// initialize
	void initialize();
	void removeEdge(int, int);
	void resetMatrix();
	inline bool isEdge(int u, int v) const { return *(matrix_Gt[u] + v); }
	bool addEdgeAndDegree(int u, int v);
	bool addEdgeNoDegree(int u, int v);
	//void addEdge(Edge);
	void addEdgeAndDegree(const string & from, const string& to,
		double cost);
	void addVertex(const string &);
	// generate the task assignment graph G_T.
	// time threshold, since we do not
	// want the task assignment process falls into 
	// an endless loop.
	bool taskAssign(const double & thre );
	void showGraph();
	void printMatrix();
	void printMatrix(std::ofstream & fs);
	void getEdges();
	void getEdges(int ** m);

	void find_4_edge_path(int s);

	void find_all_4_edge_paths(bool isDisplay, std::ofstream & fs);
	void find_all_4_edge_cycles(bool isDisplay, std::ofstream & fs);
	// with backtrack mechanism;
	void runDFS(int u, bool isDisplay, std::ofstream & fs);
	void runDFS4Cycle(int u, int des, bool isDisplay, std::ofstream & fs);
	double ** returnG_T();
	void myGraph::copyG_T(bool **);
	vector<vector<double>> myGraph::returnG_T3();
};


inline int rand_max_min_include(int min, int max){
	//srand (time(NULL));
	return rand() % (max - min + 1) + min;
}


typedef std::map<int, Edge> PATH_1_;
typedef std::vector<int> PATH_2_;

inline void printPath(int pathID, PATH_1_ & p, std::ofstream & fs){
	std::cout << "Path " << pathID << ": ";
	fs << "Path " << pathID << ": ";
	for (PATH_1_::iterator i = p.begin(); i != p.end(); ++i){
		std::cout << i->second.u << " --- " << i->second.v << ", ";
		fs << i->second.u << " --- " << i->second.v << ", ";
	}
}

inline void printPath(int pathID, PATH_2_ & p, std::ofstream & fs){
	std::cout << "Path " << pathID << ": ";
	fs << "Path " << pathID << ": ";
	for (PATH_2_::iterator i = p.begin(); i != p.end() - 1; ++i){
		std::cout << *i << " --- ";
		fs << *i << " --- ";
	}

	for (PATH_2_::iterator i = p.end() - 1; i != p.end(); ++i){
		std::cout << *i << "\n";
		fs << *i << "\n";
	}
}



#endif