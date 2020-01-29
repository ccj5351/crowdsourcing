/**
* @file: main_GT_Generator.cpp
* @brief: task assignment graph GTGenerator.
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

/* Usage:
 * Input:
 * argv[1]: vertexNum, e.g., = "200"
 * argv[2]: picking_ratio, e.g., = "0.5"
 * argv[3]; file_path. e.g., = "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/"
 */

#define _GT_Generator_CPP_ /*comment this line, or*/
#ifndef _GT_Generator_CPP_ /*change ifndef to ifdef, to make this main function work.*/
#define _CRT_SECURE_NO_WARNINGS

#include "myGraph.hpp"
#include "readCSVFiles.hpp"
#include <ctime>
#include "makeDirectory.hpp"
#include "TSPalgorithm.hpp"
#include "matrix_generator.hpp"
#include "transition.hpp"
#include <chrono> /*time*/
#include <fstream>
#include <iomanip> /*IO Manipulators*/
#include "task_assign.hpp"
#include <cstdio> /*printf*/
#include "Experiment.hpp"

using namespace std;


int main(int argc, char * argv[]){

	string s_vertexNum = argv[1]; // argv[3];
	string s_ratio = argv[2]; // argv[4];
	string file_path = argv[3];

	

	double ratio = stod(s_ratio);
	int vertexNum = stoi(s_vertexNum);
	int d = 0;
	if (ratio >= 1.0)
		d = vertexNum - 1;
	if (0 < ratio && ratio < 1.0)
		d = (int)((double)vertexNum * ratio);

	
	cout << "Generate GT : v = " << vertexNum << ", d = " << d
		<< ", ratio = " << ratio << endl;
	//generate GT
	bool ** g_t = new bool*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		g_t[i] = new bool[vertexNum];
	}

	set_GT_False(vertexNum, g_t);
	// if not successful;
	if (!task_assignment_graph(g_t, vertexNum, d)){
		return -1;
	}

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string(now->tm_hour)
		+ "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec);

	while (s_vertexNum.length() < 5){
		s_vertexNum = "0" + s_vertexNum;
	}

	string tem_name = file_path + "GT-v" + s_vertexNum + "-d" + to_string(d) + "-r-"
		+ s_ratio + "-" + cur_time + ".txt";

	std::ofstream of1(tem_name.c_str(), std::ofstream::out);

	if (!of1.is_open())
		cout << "Not Open file " << tem_name << ".\n";

	// save the task assignment graph.
	of1 << vertexNum << "\t" << ratio << "\t" << d << "\n";
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum - 1; ++j){
			of1 << g_t[i][j] << "\t";
		}
		of1 << g_t[i][vertexNum - 1] << "\n";
	}

	of1.close();
	// release memory;
	delete_2D_Array<bool>(vertexNum, g_t);
	return 0;
}
#endif