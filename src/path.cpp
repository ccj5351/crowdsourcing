/**
* @file: path.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "path.hpp"
void Path::vertex2Path(int * vertices){
		// the iterator constructor can also be used to construct from arrays:
		v_nodes = vector<int>(vertices, vertices + vertexNum);
	}


void Path::printPath(ullong & pathID){
		printf(" >- Path %llu : ", pathID);
		//std::cout << " >- Path " << pathID << " : ";
		for (uint i = 0; i < vertexNum - 1; i++)
			std::cout << "v" << v_nodes[i] << " --> ";
		std::cout << "v" << v_nodes[vertexNum - 1] << ", score = " << score << "\n";
	}

void Path::printPath(std::ofstream  & of){
	of << " >- Path : ";
	of.precision(12);
	if (v_nodes.size() > 0){
	for (uint i = 0; i < vertexNum - 1; i++)
		of << v_nodes[i] << ", ";
	of << v_nodes[vertexNum - 1] << ", score = " << score << "\n";
	printf("score = %.10f\n", score);
	}
	else{
		of << "No Hamiltonian path exists!\n";
	}
}

/* No assignment (=) operator is defined for the class.
* If you define any assignment operator that takes the
* class as a parameter, the compiler cannot generate a
* default assignment operator for that class.
* Assignment operators are not inherited by derived classes.
* You must explicitly define an assignment operator
* for each class.
*/
Path & Path:: operator= (const Path & pSource){
	pathLen = pSource.pathLen;
	vertexNum = pSource.vertexNum;
	// do the copy
	 score = pSource.score;
	 v_nodes = std::vector < int>(pSource.v_nodes);
	// return the existing object
	return *this;
}

bool isLarger(const id_weights_pair_type & arg1, const id_weights_pair_type & arg2){
	return (arg1.second > arg2.second);
}

