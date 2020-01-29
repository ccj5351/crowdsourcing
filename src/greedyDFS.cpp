/**
* @file: greedyDFS.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "greedyDFS.hpp"

//constructor 
Greedy_DFS_Path::Greedy_DFS_Path(double * * p_g,
	const int & vNum, // vertex number;
	const bool & isDisplay_,
	const bool & isSaveData_
	) :V(vNum), isDisplay(isDisplay_),
	isSaveData(isSaveData_){
	p_Gp = p_g;
	adj = new list<std::pair<int, double>>[vNum];
	for (int i = 0; i < vNum; ++i)
		for (int j = 0; j < vNum; ++j)
			if (p_g[i][j] != 0) // non-zero weights means edges;
				// Add j to i's list.
				adj[i].push_back(std::make_pair(j, p_g[i][j]));
	v_Colors = vector<VertexColors>(vNum, WHITE);
}

//destructor
Greedy_DFS_Path::~Greedy_DFS_Path(){
	for (int i = 0; i < V; ++i){
		list<pair<int, double>>().swap(adj[i]);
	}
	delete[] adj;
	vector<VertexColors>().swap(v_Colors);
	vector<ScorePathPAIR>().swap(v_HP);
}


void Greedy_DFS_Path::sortAdjListWeights(){
	// descending-sort the weights of
	// the nodes adjacent to vertex u;
	for (int u = 0; u < V; ++u)
		adj[u].sort(desSort);
}

// set color 
void Greedy_DFS_Path::setlabel(const int & u,
	const VertexColors & color){
	v_Colors[u] = color;
}

void Greedy_DFS_Path::PathDFS_version1(const int & s, std::ofstream  & of){
	// Mark the current node 
	// as visited and print it.
	if (isDisplay){
		cout << "set vertex " << s << " GRAY.\n";
		of << "set vertex " << s << " GRAY.\n";
	}
	setlabel(s, GRAY);
	s_tempPath.push(s);
	//temp_max_path_leng = std::max(int(s_tempPath.size()), 
	//	temp_max_path_leng);

	// save the currently longest path;
	// if at the end no HP occurs,
	// at least we should save the longest path.
	if (temp_max_path_leng < s_tempPath.size()){
		temp_max_path_leng = s_tempPath.size();
		s_longestPath = std::stack<int>(s_tempPath);
	}


	if (isDisplay){
		cout << " vertex " << s
			<< " gets pushed.\n";
		of << " vertex " << s
			<< " gets pushed.\n";
	}
	// if we have got an HP path;
	if (s_tempPath.size() > V - 1){
		std::cout << "-------------\n"
			<< "-------------\n";
		of << "-------------\n"
			<< "-------------\n";
		std::cout << "Already an HP.\n";
		of << "Already an HP.\n";
		cout << s_tempPath.top() << endl;
		of << s_tempPath.top() << endl;
		std::cout << "save the HP.\n";
		of << "save the HP.\n";
		printStack<int>(s_tempPath, of);
		v_HP.push_back(make_pair(getScore(), v_tempPath));
		// backtrack to the parent level
		// for another new HP.
		s_tempPath.pop();
		setlabel(s, WHITE);
		//**************
		// set flag
		//**************
		isHP = true;
		if (isDisplay){
			cout << " reset vertex " << s
				<< " WHITE, for another HP.\n";
		}
		return;
	}

	// Recur for all the vertices adjacent to this vertex
	for (auto i = adj[s].begin(); i != adj[s].end(); ++i){
		// if unexplored vertex;
		int next = (*i).first;
		if (v_Colors[next] == WHITE){
			if (isDisplay){
				std::cout << s << " ----> " << next << endl;
				of << s << " ----> " << next << endl;
			}
			//*********************
			// add something
			//*********************
			if (isHP)
				return;
			else
				PathDFS_version1(next, of);
			//*********************
			// add something
			//*********************
		}
		else{
			if (isDisplay){
				// vertex next gets popped.
				std::cout << s << " ----> " << next <<
					" -- backtracking to --> "
					<< s << ", due to " << next <<
					" has been visited before.\n";
				// vertex next gets popped.
				of << s << " ----> " << next <<
					" -- backtracking to --> "
					<< s << ", due to " << next <<
					" has been visited before.\n";
			}
		}
	}
	if (isDisplay){
		cout << "set vertex " << s << " BLACK.\n";
		of << "set vertex " << s << " BLACK.\n";
	}
	setlabel(s, BLACK);
	//vertex s gets popped;
	s_tempPath.pop();
	if (isDisplay) {
		cout << " vertex " << s
			<< " gets popped.\n";
		of << " vertex " << s
			<< " gets popped.\n";
	}
	if (s_tempPath.size() < V){
		setlabel(s, WHITE);
		if (isDisplay){
			cout << " reset vertex " << s
				<< " WHITE, due to no-HP.\n";
			of << " reset vertex " << s
				<< " WHITE, due to no-HP.\n";
		}
	}

}

void Greedy_DFS_Path::PathDFS_version2(const int & s, std::ofstream  & of){
	// Mark the current node 
	// as visited and print it.
	if (isDisplay){
		cout << "set vertex " << s << " GRAY.\n";
		of << "set vertex " << s << " GRAY.\n";
	}
	setlabel(s, GRAY);
	s_tempPath.push(s);
	// save the currently longest path;
	// if at the end no HP occurs,
	// at least we should save the longest path.
	if (temp_max_path_leng < s_tempPath.size()){
		temp_max_path_leng = s_tempPath.size();
		s_longestPath = std::stack<int>(s_tempPath);
	}
	if (isDisplay){
		cout << " vertex " << s
			<< " gets pushed.\n";
		of << " vertex " << s
			<< " gets pushed.\n";
	}
	// if we have got an HP path;
	//old condition: if (s_tempPath.size() > V - 1){
	if (s_tempPath.size() > V - 1 || isHP){
		std::cout << "-------------\n"
			<< "-------------\n";
		of << "-------------\n"
			<< "-------------\n";
		std::cout << "Already an HP. Save it.\n";
		of << "Already an HP. Save it.\n";
		//cout << s_tempPath.top() << endl;
		//of << s_tempPath.top() << endl;
		
		// save the path in stack s_tempPath
		// to vector v_tempPath;
		savePrintStack<int>(s_tempPath, of);
		v_HP.push_back(make_pair(getScore(), v_tempPath));
		//**************
		// set flag
		//**************
		isHP = true;

		//*********************
		// add something
		//*********************
		// if need exhaustive search
		// then uncomment the following.
		// backtrack to the parent level
		// for another new HP.
		// /* 
		//  * s_tempPath.pop();
		//  * setlabel(s, WHITE);
		//  * if (isDisplay){
		//	*   cout << " reset vertex " << s
		//	*	 << " WHITE, for another HP.\n";
		//  * }
		return;
	} /*end of Already an HP*/

	// if an HP is found, return;
	// else, Recursive for all the vertices adjacent to 
	// this vertex;
	if (isHP)
		return;
	// Recursive for all the vertices adjacent to 
	// this vertex;
	for (auto i = adj[s].begin(); i != adj[s].end(); ++i){
		// if unexplored vertex;
		int next = (*i).first;
		if (v_Colors[next] == WHITE){
			if (isDisplay){
				std::cout << s << " ----> " << next << endl;
				of << s << " ----> " << next << endl;
			}
			//*********************
			// add something
			//*********************
			// if need exhaustive search
			// then comment the if ... else ...
			// directly use: PathDFS(next);
			if (isHP)
				return;
			else
				PathDFS_version2(next, of);
			//*********************
			// add something
			//*********************
		}
		else{
			if (isDisplay){
				// vertex next gets popped.
				std::cout << s << " ----> " << next <<
					" -- backtracking to --> "
					<< s << ", due to " << next <<
					" has been visited before.\n";
				// vertex next gets popped.
				of << s << " ----> " << next <<
					" -- backtracking to --> "
					<< s << ", due to " << next <<
					" has been visited before.\n";
			}
		}
	} /*end of recursive call PathDFS()
	  for all the vertices adjacent to this vertex*/
	
	// if an HP is found, return;
	// else, Recursive for all the vertices adjacent to 
	// this vertex;
	if (isHP)
		return;

	if (isDisplay){
		cout << "set vertex " << s << " BLACK.\n";
		of << "set vertex " << s << " BLACK.\n";
	}
	setlabel(s, BLACK);
	//vertex s gets popped;
	s_tempPath.pop();
	if (isDisplay) {
		cout << " vertex " << s
			<< " gets popped.\n";
		of << " vertex " << s
			<< " gets popped.\n";
	}
	if (s_tempPath.size() < V){
		setlabel(s, WHITE);
		if (isDisplay){
		cout << " reset vertex " << s
			<< " WHITE, due to no-HP.\n";
		of << " reset vertex " << s
			<< " WHITE, due to no-HP.\n";
		}
	}

}

// search HPs beginning with each vertex;
void Greedy_DFS_Path::GreedyPathDFS(std::ofstream  & of){
	//sort the weights in the descending order.
	sortAdjListWeights();
	isHP = false;
	for (int u = 0; u < V; ++u){
		// Mark all the vertices as not visited
		for (int i = 0; i < V; i++){
			v_Colors[i] = WHITE;
		}
		isHP = false;
		cout << "*************************\n"
			<< "*************************\n"
			<< " Start search PathDSF(" << u << "):\n";

		of << "*************************\n"
			<< "*************************\n"
			<< " Start search PathDSF(" << u << "):\n";
		while (!s_tempPath.empty()){
			s_tempPath.pop();
		}
		temp_max_path_leng = 0;
		// Call the recursive function 
		// to print DFS traversal.
		PathDFS_version2(u, of);
		cout << "*************************\n"
			<< "*************************\n"
			<< " Finish GreedyPathDFS(" << u
			<< "), finding a path of length "
			<< temp_max_path_leng << ", i.e.,\n";

		of << "*************************\n"
			<< "*************************\n"
			<< " Finish GreedyPathDFS(" << u
			<< "), finding a path of length "
			<< temp_max_path_leng << ", i.e.,\n";
		// print the possibly longest path;
		// maybe this path is an HP,
		// also it might be a non-HP.
		printStack<int>(s_longestPath, of);

	}
	if (isDisplay){
		cout << "*************************\n"
			<< "*************************\n"
			<< "Scanning all the possible paths.\n";
		
		of  << "*************************\n"
			<< "*************************\n"
			<< "Scanning all the possible paths.\n";
	}

}

// search an HP beginning with some specific vertex;
void Greedy_DFS_Path::GreedyPathDFS(const int & v_idx, std::ofstream  & of){
	//sort the weights in the descending order.
	sortAdjListWeights();
	isHP = false;
	for (int u = v_idx; u < v_idx + 1; ++u){

		// Mark all the vertices as not visited
		for (int i = 0; i < V; i++){
			v_Colors[i] = WHITE;
		}
		isHP = false;
		cout << "*************************\n"
			<< "*************************\n"
			<< " Start GreedyPathDFS(" << u << "):\n";

		of << "*************************\n"
			<< "*************************\n"
			<< " Start GreedyPathDFS(" << u << "):\n";
		while (!s_tempPath.empty()){
			s_tempPath.pop();
		}
		temp_max_path_leng = 0;
		// Call the recursive function 
		// to print DFS traversal.
		cout << " Calling function PathDFS_version2(..)\n";
		of << " Calling function PathDFS_version2(..)\n";
		PathDFS_version2(u, of);

		//cout << " Calling function PathDFS_version1(..)\n";
		//of << " Calling function PathDFS_version1(..)\n";
		//PathDFS_version1(u);
		
		cout << "*************************\n"
			<< "*************************\n"
			<< " Finish GreedyPathDFS(" << u 
			<< "), finding a path of length "
			<< temp_max_path_leng << ", i.e.,\n";

		of << "*************************\n"
			<< "*************************\n"
			<< " Finish GreedyPathDFS(" << u
			<< "), finding a path of length "
			<< temp_max_path_leng << ", i.e.,\n";
		// print the possibly longest path;
		// maybe this path is an HP,
		// also it might be a non-HP.
		printStack<int>(s_longestPath, of);
	}

}

// descending sorting
bool desSort(const pair<int, double> & p1,
	const pair<int, double> & p2){
	return (p1.second > p2.second);
}

template < typename T >
void Greedy_DFS_Path::savePrintStack(const std::stack<T>& stk, std::ofstream  & of){
	// colon ":" :
	// This is the normal syntax for inheritance of classes/structs, 
	// here cheat is inherited from std::stack <T>.
	struct cheat : std::stack <T> { using std::stack<T>::c; };
	if (!v_tempPath.empty())
		v_tempPath.clear();
	const auto& seq = static_cast<const cheat&>(stk).c;
	for (const auto& v : seq) {
		v_tempPath.push_back(v);
		std::cout << v << ' ';
		of << v << ' ';
	}
	std::cout << '\n';
	of << '\n';
}

template < typename T >
void Greedy_DFS_Path::printStack(const std::stack<T>& stk, std::ofstream  & of){
	// colon ":" :
	// This is the normal syntax for inheritance of classes/str
	struct cheat : std::stack < T > { using std::stack<T>::c; };
	const auto& seq = static_cast<const cheat&>(stk).c;
	for (const auto& v : seq) {
		std::cout << v << ' ';
		of << v << ' ';
	}
	std::cout << '\n';
	of << '\n';
}

double Greedy_DFS_Path::getScore(){
	double score = 1.0;
	int s = v_tempPath.size();
	for (int i = 0; i < s-1; ++i){
		score *= p_Gp[v_tempPath[i]][v_tempPath[i + 1]];
	}
	return score;
}

void Greedy_DFS_Path::printHPs(std::ofstream &  of){
	if (!v_HP.empty()){
		cout << " HPs include the following:\n";
		of << " HPs include the following:\n";
		int p_idx = 1;
		for (auto i = v_HP.begin(); i != v_HP.end();
			++i, p_idx++){
			cout << "Hamiltonian Path " << p_idx << " :\n"
				<< " score = " << i->first << ", ";
			of << "Hamiltonian Path " << p_idx << " :\n"
				<< " score = " << i->first << ", ";
			auto j = i->second.begin();
			for (; j != i->second.end() - 1; ++j){
				if (isDisplay) cout << *j << " --> ";
				if (isSaveData) of << *j << " --> ";
			}
			if(isDisplay)
				cout << *j<< endl;
			if (isSaveData)
				of  << *j << endl;
		}
	}
	else{
		cout << "No HP exist!\n";
		of << "No HP exist!\n";
	}
}