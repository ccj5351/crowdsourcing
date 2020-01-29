/**
* @file: graph_ta.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "graph_ta.hpp"
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>
#include <utility>
#include "MyUtility.hpp"
#include "readCSVFiles.hpp"
#include <chrono>

using namespace std;

#define NON_NODE -1


Graph::Graph(const string& path) {
	map<int, double> m;

	this->M = readCSVFile(path);
	this->n = (int) this->M.size();
	for (int i = 0; i < this->n; i++) {
		m = map<int, double>();
		for (int j = 0; j < this->n; j++) {
			if (M[i][j] != 0)
				m[j] = M[i][j];
		}
		this->Adj.push_back(m);
	}
}

Graph::Graph(const std::vector<std::vector<double> > & gp,
	const int & N) {
	this->M = gp;
	this->n = N;
	map<int, double> m;
	for (int i = 0; i < this->n; i++) {
		m = map<int, double>();
		for (int j = 0; j < this->n; j++) {
			if (M[i][j] != 0)
				m[j] = M[i][j];
		}
		this->Adj.push_back(m);
	}
}


Graph::Graph(const string& path, const int& N) {
	map<int, double> m;

	this->M = readCSVFile(path);

	for (int i = 0; i < N; i++) {
		M[i].erase(M[i].begin() + N, M[i].end());
		M[i].shrink_to_fit();
	}
	M.erase(M.begin() + N, M.end());
	M.shrink_to_fit();
	this->n = N;

	for (int i = 0; i < this->n; i++) {
		m = map<int, double>();
		for (int j = 0; j < this->n; j++) {
			if (M[i][j] != 0)
				m[j] = M[i][j];
		}
		this->Adj.push_back(m);
	}

	cout << "n: " << this->n << endl;
	cout << "M: " << this->M.size() << endl;
	cout << this->M[0].size() << endl;
	cout << "Adj: " << this->Adj.size() << endl;
}


void
DFS_Visit(const int& previous_node, const int& current_node, vector<int>& history, double& probability, const std::vector<int>& nodes, const vector<map<int, double> >& Adj, const vector<vector<double> >& M, map<vector<int>, double>& permutations) {
	map<int, double> adj;
	int next_node;
	double transition_probability;

	//add start_node to history
	history.push_back(current_node);
	if (previous_node != NON_NODE)
		probability *= M[previous_node][current_node];

	//if history.size = nodes.size
	if (history.size() == nodes.size()) {
		//add history to permutations
		permutations[history] = probability;
	}
	else {
		//foreach adj node of the start_node that is in nodes but not history
		adj = Adj[current_node];
		for (map<int, double>::iterator it = adj.begin(); it != adj.end(); it++) {
			//DFS_Visit
			next_node = it->first;
			transition_probability = it->second;
			if ((std::find(nodes.begin(), nodes.end(), next_node) != nodes.end()) &&
				(std::find(history.begin(), history.end(), next_node) == history.end()))
				DFS_Visit(current_node, next_node, history, probability, nodes, Adj, M, permutations);
		}
	}

	//history.pop
	history.pop_back();
	if (previous_node != NON_NODE)
		probability /= M[previous_node][current_node];
}


std::map<std::vector<int>, double>
Graph::generate_permutations(const std::vector<int>& nodes) {
	map<vector<int>, double> permutations;
	double probability;
	vector<int> history;

	probability = 1;
	//foreach node
	for (vector<int>::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
		//DFS_Visit(NON_NODE, node)
		DFS_Visit(NON_NODE, *it, history, probability, nodes, this->Adj, this->M, permutations);
	}

	return permutations;
}


struct classcomp {
	bool operator() (const double& lhs, const double& rhs) const
	{
		return lhs > rhs;
	}
};


std::vector<vector<int> >
Graph::TA() {

	vector<int> path, subpath;
	vector<vector<int> > result;
	map<vector<int>, double> m;
	double threshold = 0, sub_probability;
	//the unique weights
	set<double, classcomp> weights;
	//the multimap between weight and index
	multimap<double, pair<int, int> > mm;
	pair<int, int> edge;
	multimap<double, pair<int, int> >::iterator it_low, it_up;
	map<vector<int>, double> permutations, path_weights;
	int sum = 0;
	threshold = 0;

	auto t0 = std::chrono::high_resolution_clock::now();
	//build a map between weight and index
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			edge = pair<int, int>(i, j);
			mm.insert(pair<double, pair<int, int> >(M[i][j], edge));
			//cout << M[i][j] << endl;
			weights.insert(M[i][j]);
		}
	}

	//find the edges with the largest weight
	for (set<double, classcomp>::iterator it = weights.begin();
		it != weights.end(); it++) {
		double weight = *it;
		//if the largest value is smaller than threshold, return
		if (pow(weight, n - 1) < threshold)
			break;

		//foreach edge
		it_low = mm.lower_bound(weight);
		it_up = mm.upper_bound(weight);
		for (multimap<double, pair<int, int> >::iterator it1 = it_low;
			it1 != it_up; it1++) {
			edge = (*it1).second;

			//enumerate all other integers
			subpath = vector<int>();
			for (int i = 0; i < n; i++) {
				if ((i != edge.first) && (i != edge.second))
					subpath.push_back(i);
			}

			permutations = this->generate_permutations(subpath);
			//foreach enumeration
			for (map<vector<int>, double>::iterator
				it2 = permutations.begin();
				it2 != permutations.end(); it2++){

				subpath = it2->first;
				sub_probability = it2->second;
				path_weights
					= this->insert_edge(edge, subpath, sub_probability);

				for (map<vector<int>, double>::iterator
					it3 = path_weights.begin();
					it3 != path_weights.end(); it3++) {
					if (it3->second > threshold) {
						threshold = it3->second;
						result.clear();
						result.push_back(it3->first);
					}
					else {
						if ((it3->second == threshold)
							&& (threshold != 0)) {
							result.push_back(it3->first);
						}
					}
				}
			}
		}

		sum += std::distance(it_low, it_up);
		cout << "finish " << sum << "/" << n*n << endl;
		cout << "threshold = " << threshold << endl;
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
	cout << "time: " << dt << endl;
	cout << "probability = " << threshold << endl;
	return result;
}

std::map<std::vector<int>, double>
Graph::insert_edge(const std::pair<int, int>& edge,
const std::vector<int>& subpath, const double& sub_probability) {
	map<vector<int>, double> paths;
	vector<int> path;
	double probability;

	//try to add edge to subpath.head
	if (this->M[edge.second][subpath[0]] != 0) {
		path = subpath;
		path.insert(path.begin(), edge.second);
		path.insert(path.begin(), edge.first);
		probability = M[edge.first][edge.second]
			* M[edge.second][subpath[0]]
			* sub_probability;
		paths[path] = probability;
	}

	//try to add edge to subpath.tail
	if (this->M[subpath[subpath.size() - 1]][edge.first] != 0) {
		path = subpath;
		path.push_back(edge.first);
		path.push_back(edge.second);
		probability = M[edge.first][edge.second]
			* M[subpath[subpath.size() - 1]][edge.first]
			* sub_probability;
		paths[path] = probability;
	}

	//try to add to mid
	for (int i = 1; i < subpath.size(); i++) {
		if ((this->M[subpath[i - 1]][edge.first] != 0)
			&& (this->M[edge.second][subpath[i]] != 0)) {
			path = subpath;
			path.insert(path.begin() + i, edge.second);
			path.insert(path.begin() + i, edge.first);
			probability = sub_probability
				/ M[subpath[i - 1]][subpath[i]]
				* M[edge.first][edge.second]
				* M[subpath[i - 1]][edge.first]
				* M[edge.second][subpath[i]];
			paths[path] = probability;
		}
	}

	return paths;
}