/**
* @file: task_assign.cpp
* @brief: task graph generation.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "task_assign.hpp"
#include <fstream>
using namespace std;

void deletElement(vector<int> &SetV, int element){
	std::vector<int>::iterator position = std::find(SetV.begin(), SetV.end(), element);
	if (position != SetV.end())
		SetV.erase(position);
}

bool deletEdge(vector<EDGE> &SetV, EDGE element){
	std::vector<EDGE>::iterator position = std::find(SetV.begin(), SetV.end(), element);
	if (position != SetV.end()){
		SetV.erase(position);
		return true;
	}
	else
		return false;
}

void printMatrix(const int & vertexNum, bool ** matrix_Gt){
	cout << endl;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			if (matrix_Gt[i][j])
				cout << 1;
			else
				cout << 0;
			if (vertexNum - 1 == j) cout << ";\n";
			else cout << ", ";
		}
	}
}

bool isDegreeSuccess(const int & vertexNum, bool ** matrix_Gt, const int & deg){
	for (int i = 0; i < vertexNum; ++i){
		int d = 0;
		for (int j = 0; j < vertexNum; ++j){
			if (matrix_Gt[i][j])
				d++;
		}
		if (d != deg){
			cout << "vertex " << i << " failed, its degree = "
				<< d << ".\n";
			return false;
		}
	}
	return true;
}

bool assign_task_graph(
	const int & n, // vertex num
	const int & d, // degree
	bool ** g_t
	){

	if (n*d % 2 != 0){
		std::cout << "vertex / degree = " << n << "/" << d
			<< " does not match, one of which should be even.\n";
		return false;
	}
	// complete graph, each node will connect with other nodes,
	if (d == n - 1){
		for (int i = 0; i < n; ++i)
			for (int j = i + 1; j < n; ++j)
			{
				g_t[i][j] = true;
				g_t[j][i] = true;
			}
		return true;
	}

	vector<int> v_node;
	for (int i = 0; i < n; ++i)
		v_node.push_back(i);

	// set zeros
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			g_t[i][j] = false;

	vector<int> v_degree(n, 0);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);
	std::shuffle(v_node.begin(), v_node.end(), engine);

	// generate this path, that is, to add edges connecting those vertices.
	for (int i = 0; i < n - 1; ++i){
		// add an edge (u, v), and update the degrees of u and v.
		int from = v_node[i];
		int to = v_node[i + 1];
		g_t[from][to] = true;
		g_t[to][from] = true;
		++v_degree[from];
		++v_degree[to];
	}

	if (d == 2){
		// add the first vertex and the last vertex, to generate a circle.
		g_t[v_node[0]][v_node[n - 1]] = true;
		g_t[v_node[n - 1]][v_node[0]] = true;
		++v_degree[v_node[0]];
		++v_degree[v_node[n - 1]];
		return true;
	}

	vector<int> SetV(n, 0);
	for (int i = 0; i < n; ++i)
		SetV[i] = i;
	/*
	if (d > 2 && d < n - 1){
	for (int i = 0; i < n; i++){
	//std::shuffle(v_node.begin(), v_node.end(), engine);
	if (v_degree[i] >= d)
	deletElement(SetV, i);
	while (v_degree[i] < d){
	int idx = rand() % SetV.size();
	int j = SetV[idx];
	while (j == i){
	idx = rand() % SetV.size();
	j = SetV[idx];
	}
	// already have edge (i,j) and (j,i);
	if (g_t[i][j] || g_t[j][i])
	continue;
	else{ // add (i,j)
	if (v_degree[j] < d){
	g_t[i][j] = true;
	g_t[j][i] = true;
	++v_degree[i];
	++v_degree[j];
	}
	// check degree
	if (v_degree[j] >= d)
	deletElement(SetV, j);
	}
	}//end of one value of i
	}//end of each node
	}
	*/

	if (d > 2 && d < n - 1){
		for (int i = 0; i < n; i++){
			if (v_degree[i] >= d){
				deletElement(SetV, i);
				continue;
			}

			std::shuffle(v_node.begin(), v_node.end(), engine);
			for (int k = 0; k < SetV.size(); ++k){
				int j = v_node[k];
				if (j == i)
					continue;

				// already have edge (i,j) and (j,i);
				if (g_t[i][j] || g_t[j][i])
					continue;
				else{ // add (i,j)
					if (v_degree[j] < d){
						g_t[i][j] = true;
						g_t[j][i] = true;
						++v_degree[i];
						++v_degree[j];
					}
				}
			}/*end of one value of i*/
		}/*end of each node*/
	}
	return isDegreeSuccess(n, g_t, d);
}

// task assignment via deleting edges,
// that is, decrease degree from v-1 to d.
// where, d <= v-1;
bool assign_task_graph_down_2_d(
	const int & n, // vertex num
	const int & d, // degree
	bool ** g_t,
	ofstream & of
	){
	if (d < 2){
		//std::cout << "degree = " << d << " < 2\n";
		of << "degree = " << d << " < 2\n";
		return false;
	}

	if (n*d % 2 != 0){
		//std::cout << "vertex / degree = " << n << "/" << d
		//	<< " does not match, one of which should be even.\n";
		of << "vertex / degree = " << n << "/" << d
			<< " does not match, one of which should be even.\n";
		return false;
	}
	// complete graph, each node will connect with other nodes,
	if (d == n - 1){
		for (int i = 0; i < n; ++i)
			for (int j = i + 1; j < n; ++j)
			{
				g_t[i][j] = true;
				g_t[j][i] = true;
			}
		return true;
	}

	vector<int> v_node;
	for (int i = 0; i < n; ++i)
		v_node.push_back(i);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);
	std::shuffle(v_node.begin(), v_node.end(), engine);
	// zeros 
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			g_t[i][j] = false;
		}

	// generate this path, that is, to add edges connecting those vertices.
	for (int i = 0; i < n - 1; ++i){
		// add an edge (u, v), and update the degrees of u and v.
		int from = v_node[i];
		int to = v_node[i + 1];
		g_t[from][to] = true;
		g_t[to][from] = true;
	}

	if (d == 2){
		// add the first vertex and the last vertex, to generate a circle.
		g_t[v_node[0]][v_node[n - 1]] = true;
		g_t[v_node[n - 1]][v_node[0]] = true;
		return true;
	}

	// initial degree is v-1;
	map<int, int> m_degrees;
	for (int i = 0; i < n; ++i){
		m_degrees[i] = n - 1;
	}

	vector<std::pair<int, int>> v_edges;
	// complete edges container;
	for (int i = 0; i < n; ++i){
		for (int j = i + 1; j < n; ++j){
			v_edges.push_back(make_pair(i, j)); // edge (u, v), where u < v;
		}
	}

	// reserve the edges in the path;
	for (int i = 0; i < n - 1; ++i){
		// edge (u, v), where u < v;
		int u = v_node[i];
		int v = v_node[i + 1];
		EDGE e = make_pair(std::min(u, v), std::max(u, v));
		if (!deletEdge(v_edges, e)){
			cout << "Cannot delete edge (" << std::min(u, v)
				<< ", " << std::max(u, v) << ").\n";
			return false;
		}
	}

	// delete other edges until the degree = d;
	int numEdgeToBeDeleted = n*(n - 1) / 2 - n*d / 2;
	int count = 0;
	vector<EDGE> v_temp_Edges;
	vector<EDGE> v_copy_edges = v_edges;
	int IterationTimes = 10000;
	bool isSucess = false;
	while (IterationTimes > 0){
		// keep deleting edges until arriving d;
		count = 0;
		v_edges = v_copy_edges;
		shuffle(v_edges.begin(), v_edges.end(), engine);

		if (!v_temp_Edges.empty())
			v_temp_Edges.clear();
		// find some edges which will be deleted later.
		for (int i = 0, j = v_edges.size() - 1; i <= j; ++i, --j){
			// edge (u, v), where u < v;
			int u1 = std::min(v_edges[i].first, v_edges[i].second);
			int v1 = std::max(v_edges[i].first, v_edges[i].second);
			int u2 = std::min(v_edges[j].first, v_edges[j].second);
			int v2 = std::max(v_edges[j].first, v_edges[j].second);

			if (m_degrees[u1] > d && m_degrees[v1] > d){
				EDGE temp_e = make_pair(u1, v1);
				v_temp_Edges.push_back(temp_e);
				// change the degree;
				(m_degrees[u1])--;
				(m_degrees[v1])--;
				count++;
			}

			if ((i != j) && m_degrees[u2] > d && m_degrees[v2] > d){
				EDGE temp_e2 = make_pair(u2, v2);
				v_temp_Edges.push_back(temp_e2);
				// change the degree;
				(m_degrees[u2])--;
				(m_degrees[v2])--;
				count++;
			}

			if (count >= numEdgeToBeDeleted)
			{
				isSucess = true;
				break;
			}
		}/*end of for-loop*/

		// remove the edges;
		for (auto j : v_temp_Edges){
			if (!deletEdge(v_edges, j)){
				//cout << "Cannot delete edge (" << j.first
				//	<< ", " << j.second << ").\n";
				return false;
			}
		}

		// clear temporary container;
		v_temp_Edges.clear();
		IterationTimes--;
		if (count >= numEdgeToBeDeleted)
		{
			isSucess = true;
			break;
		}

	}/*end of while-loop*/

	if (isSucess){
		for (auto i : v_edges){
			int u = i.first;
			int v = i.second;
			g_t[u][v] = true;
			g_t[v][u] = true;
		}
		return isDegreeSuccess(n, g_t, d);
	}
	else
		return false;
}
