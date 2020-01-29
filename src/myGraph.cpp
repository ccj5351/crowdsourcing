/**
* @file: myGraph.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "myGraph.hpp"

using namespace std;
//using namespace Eigen;

void myGraph::addVertex(const string & name){
	vmap::iterator itr = work.begin();
	itr = work.find(name);
	if (itr == work.end()){
		vertex * v = new vertex(name);
		work[name] = v;
		return;
	}
	std::cout << "\nVertex " << name << " already exists!\n";
}

vector<vector<double>> myGraph::returnG_T3(){
	vector<vector<double>> g_t(vertexNum, vector<double>(vertexNum, .0));
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			g_t[i][j] = (matrix_Gt[i][j] == true) ? 1.0 : 0.0;
		}
	}
	return g_t;
}

double ** myGraph::returnG_T(){
	double **pp_g_t = new double*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		pp_g_t[i] = new double[vertexNum];
	}
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			pp_g_t[i][j] = (matrix_Gt[i][j] == true) ? 1.0 : 0.0;
		}
	}
	return pp_g_t;
}

void myGraph::copyG_T(bool ** pp_g_t ){
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			pp_g_t[i][j] = matrix_Gt[i][j];
		}
	}
}

bool myGraph::addEdgeNoDegree(int from, int to){
	if (from == to)
		return false;
	if (isDirected){ // directed graph;
		if (matrix_Gt[from][to]) // this edge already exists!
			return false;
		else{
			matrix_Gt[from][to] = true;
			return true;
		}
	}
	else{ // undirected graph;
		// this edge already exists!
		if (matrix_Gt[from][to] || matrix_Gt[to][from])
			return false;
		else {
			matrix_Gt[from][to] = true;
			matrix_Gt[to][from] = true;
			return true;
		}
	}
}

bool myGraph::addEdgeAndDegree(int from, int to){

	if (from == to)
		return false;
	if (isDirected){ // directed graph;
		if (matrix_Gt[from][to]) // this edge already exists!
			return false;
		else{
			matrix_Gt[from][to] = true;
			v_degree[from] += 1;
			v_degree[to] += 1;
			return true;
		}
	}
	else{ // undirected graph;
		// this edge already exists!
		if (matrix_Gt[from][to] || matrix_Gt[to][from])
			return false;
		else {
			matrix_Gt[from][to] = true;
			matrix_Gt[to][from] = true;
			v_degree[from] += 1;
			v_degree[to] += 1;
			return true;
		}
	}
}

void myGraph::removeEdge(int from, int to){
	*(matrix_Gt[from] + to) = false;
}


// generate the task assignment graph G_T.
bool myGraph::taskAssign(
	const double & thre // time threshold, since we do not
	// want the task assignment process falls into 
	// an endless loop.
	){

	double delta_time = .0; // duration in seconds;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	// 1.e-9: nano-seconds to seconds;
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

	// for each vertex, randomly add d_min edges for it.
	int rand_v = 0, rand_u = 0, idx_u = 0, idx_v = 0;
	int temp_degree = 0;

	// set vertex set V = {0, 1, ..., |V|-1}:
	vector <int> v_temp_vertex(vertexNum, 0);
	for (int i = 0; i < vertexNum; ++i)
		v_temp_vertex[i] = i;

	// firstly, randomly generate a Hamiltonian path 
	// consisting of all the vertices.
	// using built-in random generator:

	/* initialize random seed: */
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);
	std::shuffle(v_temp_vertex.begin(), v_temp_vertex.end(), engine);
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	// print the path content:
	std::cout << "The randomly permuted vertices contains: ";
	for (auto i : v_temp_vertex)
		std::cout << ' ' << i;
	std::cout << '\n';
#endif __RELEASE_MODE_CORWDSOURCING_PROJECT_
	// generate this path, that is, to add edges connecting those vertices.
	for (std::vector<int>::iterator it = v_temp_vertex.begin();
		it != v_temp_vertex.end() - 1; ++it){
		// add an edge (u, v), and update the degrees of u and v.
		addEdgeAndDegree(*it, *(it + 1));
	}


	if (this->degree_min == 1){
		std::cout << "For an HP, the degree should be >= 2.\n";
		return false;
	}

	else if (this->degree_min == 2){
		// add the first vertex and the last vertex, to generate a circle.
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
		std::cout << "degree = 2, to make a circle:\n";
#endif
		addEdgeAndDegree(v_temp_vertex[0], v_temp_vertex[vertexNum - 1]);
		return true;
	}

	else if (this->degree_min > 2 && this->degree_min < vertexNum - 1){
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
		std::cout << "2 < degree = " << degree_min << " <= " 
			<< vertexNum - 1 << "(i.e., V -1)\n";
#endif
		//int temp_size = v_temp_vertex.size();
		while (v_temp_vertex.size() > 3){
			// check the time;
			t2 = std::chrono::high_resolution_clock::now();
			dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
			if (dt > thre){
				std::cout << "Time too long, break the loop!\n";
				return false;
			}
			std::shuffle(v_temp_vertex.begin(), v_temp_vertex.end(), engine);
			int min_d_vIdx = 0;
			// find the vertex with the least degrees.
			for (int i = 0; i < v_temp_vertex.size(); ++i){
				if (v_degree[v_temp_vertex[min_d_vIdx]]
			> v_degree[v_temp_vertex[i]])
			min_d_vIdx = i;
			}

			idx_u = min_d_vIdx;
			idx_v = rand() % v_temp_vertex.size();
			while (idx_v == idx_u){
				idx_v = rand() % v_temp_vertex.size();
			}
			rand_u = v_temp_vertex[idx_u];
			rand_v = v_temp_vertex[idx_v];
			//add an edge between vertex rand_u and rand_v.
			addEdgeAndDegree(rand_u, rand_v);

			// if degree requirement of vertex rand_u is satisfied,
			// and degree requirement of vertex rand_v is satisfied.
			// then remove vertex rand_u and rand_v from the vertex set.
			if (v_degree[rand_u] >= degree_min && v_degree[rand_v] >= degree_min){
				v_temp_vertex.erase(v_temp_vertex.begin()
					+ std::min(idx_u, idx_v));
				v_temp_vertex.erase(v_temp_vertex.begin()
					+ std::max(idx_u, idx_v) - 1);
			}
			// if degree requirement of vertex rand_u is satisfied,
			// but degree requirement of vertex rand_v is not satisfied.
			// then remove vertex rand_u from the vertex set.
			if (v_degree[rand_u] >= degree_min && v_degree[rand_v] < degree_min)
				v_temp_vertex.erase(v_temp_vertex.begin() + idx_u);

			// if degree requirement of vertex rand_v is satisfied,
			// but degree requirement of vertex rand_u is not satisfied.
			// then remove vertex rand_v from the vertex set.
			if (v_degree[rand_u] < degree_min && v_degree[rand_v] >= degree_min)
				v_temp_vertex.erase(v_temp_vertex.begin() + idx_v);

		}/*end of while*/

#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
		std::cout << "The last few steps, here temp_vertex_set_size = "
			<< v_temp_vertex.size() << endl;
#endif
		while (v_temp_vertex.size() > 0){
			if (v_temp_vertex.size() == 3){
				// pick the vertex with the minimum degrees out of
				// those three left vertices.

				int min_d_vIdx = 0;
				// find the vertex with the least degrees.
				for (int i = 0; i < v_temp_vertex.size(); ++i){
					if (v_degree[v_temp_vertex[min_d_vIdx]]
				> v_degree[v_temp_vertex[i]])
				min_d_vIdx = i;
				}

				std::cout << "min_d_vertex = " << v_temp_vertex[min_d_vIdx] << endl;
				// try to add edges starting from the vertex with 
				// the least degree. 

				// sometimes, the last 3 nodes cannot be connected to each 
				// another, that means the vertex size cannot shrink,
				// thus, the while-loop will never stop, result in an
				// endless loop. Pay attention to this possible endless loop.
				bool shrinkFlag = false;
				for (auto i : v_temp_vertex){
					if (!addEdgeAndDegree(i, v_temp_vertex[min_d_vIdx]))
					{
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_						
						std::cout << ", cannot connect " << i << " and "
							<< v_temp_vertex[min_d_vIdx] << "\n";
#endif
					}
					else
						shrinkFlag = true;
				}

				int temp_size = v_temp_vertex.size();
				for (int i = 0; i < v_temp_vertex.size(); i++)
					if (v_degree[v_temp_vertex[i]] >= degree_min){
					v_temp_vertex.erase(v_temp_vertex.begin() + i);
					// due to removal of one element,
					// adjust the index by subtracting 1;
					i--;
					}
				// no edge between the last 3 vertices,
				// resulting an endless loop, we have to return false;
				if (!shrinkFlag){
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
					std::cout << "The last 3 vertices cannot generate"
					<< " any further edges. No shrink, thus resulting"
				    << " in an endless loop, so return false right now.\n";
#endif
					return false;
				}
			}
			else if (v_temp_vertex.size() == 2){
				if (addEdgeAndDegree(v_temp_vertex[0], v_temp_vertex[1])){
					// remove those two nodes, since their degree >= degree_min.
					if (v_degree[0] >= degree_min && v_degree[1] >= degree_min)
						v_temp_vertex.erase(v_temp_vertex.begin(),
						v_temp_vertex.begin() + 2);
					// remove one nodes, since its degree >= degree_min.
					if (v_degree[0] >= degree_min && v_degree[1] < degree_min)
						v_temp_vertex.erase(v_temp_vertex.begin());
					// remove one nodes, since its degree >= degree_min.
					if (v_degree[0] < degree_min && v_degree[1] >= degree_min)
						v_temp_vertex.erase(v_temp_vertex.begin() + 1);
				}
				else{
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
					std::cout << "cannot connect the left 2 vertices!\n";
#endif
					return false;
				}
			}

			else{
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
				std::cout << "The left vertices are " << v_temp_vertex.size()
					<< " vertices, that is, neither"
					<< " 3 vertices nor 2 vertices at all.\n";
#endif
				return false;
			}
		}
		if (v_temp_vertex.size() == 0)
			std::cout << " Task assignment succeed!\n";
		return true;
	}/* 2 < degree < v-1 */

	// complete graph
	// thus, connect each i and j, where i != j;
	else if (this->degree_min == vertexNum - 1)
	{
		for (int i = 0; i < vertexNum; ++i)
			for (int j = i + 1; j < vertexNum; ++j)
				addEdgeAndDegree(i, j);
		std::cout << " Complete Graph, Task assignment succeed!\n";
		return true;
	}
	else{
		std::cout << " Error! Degree d = "<< degree_min
			<< " is not correct at all!\n";
		return false;
	}
}

void myGraph::showGraph(){
	for (int i = 0; i < vertexNum; ++i){
		std::cout << "Vertex " << i << " connects with "
			<< v_degree[i] << " vertices, including: \n";
		for (int j = 0; j < vertexNum; ++j)
			if (isEdge(i, j))
				cout << j << ", ";
		cout << endl;
	}
}

// initialize
void myGraph::initialize(){
	resetMatrix();
	for (auto i = 0; i < v_degree.size(); ++i)
		v_degree[i] = 0;
}
// reset the value to false
void myGraph::resetMatrix(){
	for (int i = 0; i < vertexNum; ++i)
		for (int j = 0; j < vertexNum; ++j)
			matrix_Gt[i][j] = false;
}
void myGraph::printMatrix(std::ofstream & of1){
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			if (matrix_Gt[i][j])
				of1 << 1;
			else
				of1 << 0;
			if (vertexNum - 1 == j) of1 << "\n";
			else of1 << ", ";
		}
	}
}


void myGraph::setSink(const int & v){
	for (auto i = 0; i < vertexNum; ++i)
		// out-degree(sink) = 0;
		matrix_Gt[v][i] = false;
}
void myGraph::setSource(const int & v){
	for (auto i = 0; i < vertexNum; ++i)
		// in-degree(source) = 0;
		matrix_Gt[i][v] = false;
}


void myGraph::setIsolateVertex(const int & v){
	for (auto i = 0; i < vertexNum; ++i){
		// in-degree(isolated vertex) = 0;
		matrix_Gt[i][v] = false;
		// out-degree(isolated vertex) = 0;
		matrix_Gt[v][i] = false;
	}
}
void myGraph::printMatrix(){
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


void myGraph::getEdges(){
	for (int i = 0; i < vertexNum; ++i){
		for (int j = i + 1; j < vertexNum; ++j){
			if (*(matrix_Gt[j] + i))
				//m_undirectedEdges.push_back(Edge(i, j));
				m_undirectedEdges[i * vertexNum + j] = Edge(i, j);

		}
	}
	cout << "Task assignment Graph G_T has "
		<< m_undirectedEdges.size() << " edges.\n";
}
void myGraph::getEdges(int ** m){
	for (int i = 0; i < vertexNum; ++i)
		v_degree[i] = 0;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			*(matrix_Gt[j] + i) = *(m[j] + i);
			if (*(m[j] + i) && (j > i)) {
				//m_undirectedEdges.push_back(Edge(i, j));
				m_undirectedEdges[i * vertexNum + j] = Edge(i, j);
				v_degree[i] += 1;
				v_degree[j] += 1;
			}
		}
	}
	cout << "Task assignment Graph G_T has " << m_undirectedEdges.size() << " edges.\n";
}


void myGraph::find_all_4_edge_paths(bool isDisplay, std::ofstream & fs){
	int pathID = 0;
	if (fs.is_open()) cout << "Open!\n";
	else cout << "not open!\n";
	for (int i = 0; i < vertexNum; ++i)
		temp_vertexSet.push_back(i);
	/* initialize random seed: */
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);
	while (!m_undirectedEdges.empty()){
		if (isDisplay)
			cout << m_undirectedEdges.size() << " edges are left.\n";
		for (int i = 0; i < vertexNum; ++i)
			v_Colors[i] = WHITE;
		// using built-in random generator:
		std::shuffle(temp_vertexSet.begin(), temp_vertexSet.end(), engine);
		if (m_undirectedEdges.size() <= path_Edge_length){
			cout << path_Edge_length <<
				" edges are left, edge-assignment task done!\n";
			if (isDisplay)
				printPath(++pathID, m_undirectedEdges, fs);
			m_undirectedEdges.clear();
		}
		else{
			runDFS(temp_vertexSet[0], isDisplay, fs);
			if (isDisplay)
				printPath(++pathID, tempPath, fs);
			/*
			if (isDisplay){
			std::cout << "Path " << pathID << ": ";
			fs << "Path " << pathID << ": ";
			for (auto i = tempPath.begin(); i != tempPath.end() - 1; ++i){
			std::cout << *i << " --- ";
			fs << *i << " --- ";
			}
			// ?????
			for (auto i = tempPath.end() - 1; i != tempPath.end(); ++i){
			std::cout << *i << "\n";
			fs << *i << "\n";
			}
			}*/
			tempPath.clear();
		}
	}
	fs.close();
}

void myGraph::runDFS(int u, bool isDisplay, std::ofstream & fs){
	// firstly add the vertex u, 
	// which will be later removed if u
	// is not satisfactory, then we backtrack
	// to u's parent.
	tempPath.push_back(u);
	v_Colors[u] = GRAY;
	//int possible_add_edge = 0;
	for (int v = 0; v < vertexNum; v++){

		if (tempPath.size() < path_Vertex_Length){
			if (isEdge(u, v) && v_Colors[v] == WHITE){
				if (isDisplay)
					std::cout << "  " << u << " ----> " << v << endl;
				fs << "  " << u << " ----> " << v << endl;

				// If a 5-vertex path has already been , 
				// then return, else keep doing dfs traversal;
				if (tempPath.size() < path_Vertex_Length){
					v_parents[v] = u;
					removeEdge(u, v);
					if (!isDirected) removeEdge(v, u);
					v_degree[u] -= 1;
					v_degree[v] -= 1;
					if (v_degree[u] <= 0){
						auto it = std::find(temp_vertexSet.begin(),
							temp_vertexSet.end(), u);
						if (it != temp_vertexSet.end()){
							temp_vertexSet.erase(it);
							if (isDisplay)
								cout << "*********************\n"
								"deg(" << u << ") <= 0, delete vertex "
								<< u << "\n*********************\n";
						}
					}

					if (v_degree[v] <= 0){
						auto it = std::find(temp_vertexSet.begin(),
							temp_vertexSet.end(), v);
						if (it != temp_vertexSet.end()){
							temp_vertexSet.erase(it);
							if (isDisplay)
								cout << "*********************\n"
								"deg(" << v << ") <= 0, delete vertex "
								<< v << "\n*********************\n";
						}
					}

					// map , erasing by key;
					// erasing by iterator is also available.
					// here we use max(u, v), and min(u, v),
					// because when saving the undirected edges,
					// we only consider the edge (i, j) or (j, i)
					// as (min(i, j), max(i, j)) = (u, v), i.e., u < v;
					if (isDisplay)
						cout << " m_undirectedEdges set erase edge "
						<< u << " --- " << v << endl;
					//possible_add_edge = v;
					m_undirectedEdges.erase(std::min(u, v)* vertexNum
						+ std::max(u, v));
					runDFS(v, isDisplay, fs);
				}
				else {
					std::cout << " cannot adding vertex " << v
						<< ", already generate a " << path_Vertex_Length
						<< "-vertex path, thus stop the search.\n";
					fs << " cannot adding vertex " << v
						<< ", already generate a " << path_Vertex_Length
						<< "-vertex path, thus stop the search.\n";
					return;
				}
			}

			else if (isEdge(u, v) && v_Colors[v] != WHITE){
				if (isDisplay){
					std::cout << "  " << u << " ----> " << v <<
						" -- backtracking to --> "
						<< u << ", due to " << v <<
						" has already been visited before." << endl;
					fs << "  " << u << " ----> " << v <<
						" -- backtracking to --> "
						<< u << ", due to " << v <<
						" has already been visited before." << endl;
				}
			}
		}
		else { // path finished!
			if (isDisplay)
				cout << " Already get a " << path_Vertex_Length
				<< "-vertex path!"
				<< " No adding edge (" << u << " , " << v
				<< "). Immediately next vertex!\n";
			if (isDisplay)
				fs << " Already get a " << path_Vertex_Length
				<< "-vertex path!"
				<< " No adding edge (" << u << " , " << v
				<< "). Immediately next vertex!\n";
			return;
		}
	}
	v_Colors[u] = BLACK;
	// undo edge delete
	if (tempPath.size() > 1 &&
		tempPath.size() < path_Vertex_Length){
		if (isDisplay)
			cout << " undo edge deletion m_undirectedEdges set re-add edge "
			<< v_parents[u] << " --- " << u << endl;
		m_undirectedEdges[std::min(u, v_parents[u])* vertexNum
			+ std::max(u, v_parents[u])] =
			Edge(std::min(u, v_parents[u]), std::max(u, v_parents[u]));

		if (addEdgeNoDegree(v_parents[u], u) && isDisplay)
			cout << " deg(" << u << ") + 1, and "
			<< " deg(" << v_parents[u] << ") + 1\n";
		else
			cout << " Error occurs for deg(" << u << ") + 1, and "
			<< " deg(" << v_parents[u] << ") + 1\n";
		// undo edge deletion,
		// that is, add 1 to the degree of the 
		// involved vertices.
		v_degree[u] += 1;
		v_degree[v_parents[u]] += 1;
		// re-add the vertex 
		for (int i = 0; i < vertexNum; ++i){
			if (v_degree[i] > 0){
				auto it = std::find(temp_vertexSet.begin(),
					temp_vertexSet.end(), i);
				if (it == temp_vertexSet.end()){
					temp_vertexSet.push_back(i);
					if (isDisplay)
						cout << "*********************\n"
						"deg(" << i << ") > 0, re-add vertex "
						<< i << "\n*********************\n";
				}
			}
		}


		// although the vertex u was added to 
		// the temporary path during some
		// previous procedure, it has to be
		// removed since u is not satisfactory, 
		// and we have to backtrack to its parent.
		auto it = std::find(tempPath.begin(), tempPath.end(), u);
		if (it != tempPath.end())
			tempPath.erase(it);

		if (isDisplay)
			cout << " undo path-vertex adding, remove vertex "
			<< u << " from tempPath.\n";
	}
	else if (1 == tempPath.size()){
		if (isDisplay)
			cout << "  tempPath.size = 1\n";
	}

	if (isDisplay)
		cout << "  Set vertex " << u << " be BLACK or VISITED."
		<< " Backtrack to its parent " << v_parents[u]
		<< " to next round DFS.\n";
	if (isDisplay)
		fs << "  Set vertex " << u << " be BLACK or VISITED."
		<< " Backtrack to its parent " << v_parents[u]
		<< " for next round DFS.\n";
}

void myGraph::find_all_4_edge_cycles(bool isDisplay, std::ofstream & fs){
	int pathID = 0;
	if (fs.is_open()) cout << "Open!\n";
	else cout << "not open!\n";
	for (int i = 0; i < vertexNum; ++i)
		temp_vertexSet.push_back(i);
	/* initialize random seed: */
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);
	while (!m_undirectedEdges.empty()){
		if (isDisplay)
			cout << m_undirectedEdges.size() << " edges are left.\n";
		for (int i = 0; i < vertexNum; ++i)
			v_Colors[i] = WHITE;
		// using built-in random generator:
		std::shuffle(temp_vertexSet.begin(), temp_vertexSet.end(), engine);
		if (m_undirectedEdges.size() <= path_Edge_length){
			if (isDisplay){
				cout << path_Edge_length <<
					" edges are left, edge-assignment task done!\n";
				printPath(++pathID, m_undirectedEdges, fs);
			}
			m_undirectedEdges.clear();
		}
		else{
			runDFS4Cycle(temp_vertexSet[0], temp_vertexSet[0], isDisplay, fs);
			if (isDisplay)
				printPath(++pathID, tempPath, fs);
			/*
			if (isDisplay){
			std::cout << "Path " << pathID << ": ";
			fs << "Path " << pathID << ": ";
			for (auto i = tempPath.begin(); i != tempPath.end() - 1; ++i){
			std::cout << *i << " --- ";
			fs << *i << " --- ";
			}
			// ?????
			for (auto i = tempPath.end() - 1; i != tempPath.end(); ++i){
			std::cout << *i << "\n";
			fs << *i << "\n";
			}
			}*/
			tempPath.clear();
		}
	}
	fs.close();
}

void myGraph::runDFS4Cycle(int u, int des, bool isDisplay, std::ofstream & fs){
	// firstly add the vertex u, 
	// which will be later removed if u
	// is not satisfactory, then we backtrack
	// to vertex u's parent.
	auto it = std::find(tempPath.begin(), tempPath.end(), u);
	if (it == tempPath.end()){
		tempPath.push_back(u);
	}
	v_Colors[u] = GRAY;
	if (u == des && tempPath.size() == cycle_Vertex_Length)
		return;

	//int possible_add_edge = 0;
	for (int v = 0; v < vertexNum; v++){

		// still use path_Vertex_Length,
		// instead of cycle_Vertex_Length.
		if (tempPath.size() < cycle_Vertex_Length){
			if (isEdge(u, v) && v_Colors[v] == WHITE){
				if (isDisplay)
					std::cout << "  " << u << " ----> " << v << endl;
				fs << "  " << u << " ----> " << v << endl;

				// If a 5-vertex path has already been , 
				// then return, else keep doing dfs traversal;
				if (tempPath.size() < cycle_Vertex_Length){
					v_parents[v] = u;
					removeEdge(u, v);
					if (!isDirected) removeEdge(v, u);
					v_degree[u] -= 1;
					v_degree[v] -= 1;
					if (v_degree[u] <= 0){
						auto it = std::find(temp_vertexSet.begin(),
							temp_vertexSet.end(), u);
						if (it != temp_vertexSet.end()){
							temp_vertexSet.erase(it);
							if (isDisplay)
								cout << "*********************\n"
								"deg(" << u << ") <= 0, delete vertex "
								<< u << "\n*********************\n";
						}
					}

					if (v_degree[v] <= 0){
						auto it = std::find(temp_vertexSet.begin(),
							temp_vertexSet.end(), v);
						if (it != temp_vertexSet.end()){
							temp_vertexSet.erase(it);
							if (isDisplay)
								cout << "*********************\n"
								"deg(" << v << ") <= 0, delete vertex "
								<< v << "\n*********************\n";
						}
					}

					// map , erasing by key;
					// erasing by iterator is also available.
					// here we use max(u, v), and min(u, v),
					// because when saving the undirected edges,
					// we only consider the edge (i, j) or (j, i)
					// as (min(i, j), max(i, j)) = (u, v), i.e., u < v;
					if (isDisplay)
						cout << " m_undirectedEdges set erase edge "
						<< u << " --- " << v << endl;
					//possible_add_edge = v;
					m_undirectedEdges.erase(std::min(u, v)* vertexNum
						+ std::max(u, v));
					runDFS4Cycle(v, des, isDisplay, fs);
				}
				else {
					std::cout << " cannot adding vertex " << v
						<< ", already generate a " << path_Vertex_Length
						<< "-vertex path, thus stop the search.\n";
					fs << " cannot adding vertex " << v
						<< ", already generate a " << path_Vertex_Length
						<< "-vertex path, thus stop the search.\n";
					return;
				}
			}

			// GRAY or Black color, then back-track
			else if (isEdge(u, v) && v_Colors[v] != WHITE){
				if (isDisplay){
					std::cout << "  " << u << " ----> " << v <<
						" -- backtracking to --> "
						<< u << ", due to " << v <<
						" has already been visited before." << endl;
					fs << "  " << u << " ----> " << v <<
						" -- backtracking to --> "
						<< u << ", due to " << v <<
						" has already been visited before." << endl;
				}
			}
		}
		else { // path finished!
			if (isDisplay)
				cout << " Already get a " << path_Vertex_Length
				<< "-vertex path!"
				<< " No adding edge (" << u << " , " << v
				<< "). Immediately next vertex!\n";
			if (isDisplay)
				fs << " Already get a " << path_Vertex_Length
				<< "-vertex path!"
				<< " No adding edge (" << u << " , " << v
				<< "). Immediately next vertex!\n";
			return;
		}
	}
	v_Colors[u] = BLACK;
	// undo edge delete
	if (tempPath.size() > 1 &&
		tempPath.size() < path_Vertex_Length){
		if (isDisplay)
			cout << " undo edge deletion m_undirectedEdges set re-add edge "
			<< v_parents[u] << " --- " << u << endl;
		m_undirectedEdges[std::min(u, v_parents[u])* vertexNum
			+ std::max(u, v_parents[u])] =
			Edge(std::min(u, v_parents[u]), std::max(u, v_parents[u]));

		if (addEdgeNoDegree(v_parents[u], u) && isDisplay)
			cout << " deg(" << u << ") + 1, and "
			<< " deg(" << v_parents[u] << ") + 1\n";
		else
			cout << " Error occurs for deg(" << u << ") + 1, and "
			<< " deg(" << v_parents[u] << ") + 1\n";
		// undo edge deletion,
		// that is, add 1 to the degree of the 
		// involved vertices.
		v_degree[u] += 1;
		v_degree[v_parents[u]] += 1;
		// re-add the vertex 
		for (int i = 0; i < vertexNum; ++i){
			if (v_degree[i] > 0){
				auto it = std::find(temp_vertexSet.begin(),
					temp_vertexSet.end(), i);
				if (it == temp_vertexSet.end()){
					temp_vertexSet.push_back(i);
					if (isDisplay)
						cout << "*********************\n"
						"deg(" << i << ") > 0, re-add vertex "
						<< i << "\n*********************\n";
				}
			}
		}


		// although the vertex u was added to 
		// the temporary path during some
		// previous procedure, it has to be
		// removed since u is not satisfactory, 
		// and we have to backtrack to its parent.
		auto it = std::find(tempPath.begin(), tempPath.end(), u);
		if (it != tempPath.end())
			tempPath.erase(it);

		if (isDisplay)
			cout << " undo path-vertex adding, remove vertex "
			<< u << " from tempPath.\n";
	}
	else if (1 == tempPath.size()){
		if (isDisplay)
			cout << "  tempPath.size = 1\n";
	}

	if (isDisplay)
		cout << "  Set vertex " << u << " be BLACK or VISITED."
		<< " Backtrack to its parent " << v_parents[u]
		<< " to next round DFS.\n";
	if (isDisplay)
		fs << "  Set vertex " << u << " be BLACK or VISITED."
		<< " Backtrack to its parent " << v_parents[u]
		<< " for next round DFS.\n";
}