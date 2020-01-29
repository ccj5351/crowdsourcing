/**
* @file: graph_ta.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef graph_hpp
#define graph_hpp

#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <utility>

class Graph {
public:
	Graph(const std::string& path);
	Graph::Graph(const std::vector<std::vector<double> > & gp, const int& N);
	//Graph::Graph(double **  gp, const int& N);
	Graph(const std::string& path, const int& N);
	std::map<std::vector<int>, double> generate_permutations(const std::vector<int>& nodes);   //generate the permutations using the graph structure
	std::map<std::vector<int>, double> insert_edge(const std::pair<int, int>& edge, const std::vector<int>& subpath, const double& sub_probability);
	std::vector<std::vector<int> > TA();

	std::vector<std::map<int, double> > Adj;//the adjacency list
	std::vector<std::vector<double> > M;
	int n;  //# of objects
};

#endif /* graph_hpp */
