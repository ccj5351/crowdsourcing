/**
* @file: main_GT_VM_GP_Generator.cpp
* @brief: task assignment graph GTGenerator.
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/


/* Usage:
 ***********
 * Input:***
 ***********
 * argv[1]: vertexNum, e.g., = "200"
 * argv[2]: picking_ratio, e.g., = "0.5"
 * argv[3]: Gauss_mean. e.g., = "0.0";// Gaussian variance to control the worker's quality.
 * argv[4]: Gauss_stddev, e.g., = "0.001";// Gaussian variance to control the worker's quality.
 * argv[5]: distribuion_type, e.g., = 1(Gaussian Distribtution); = 0(Uniform distribution)
 * argv[6]: workerNum, e.g., = "50"; the number of workers.
 * argv[7]: output_file_path, e.g., = "E:/CrowdSourcing/ICDE-2017/"
 
 ************
 *Output: a txt file, including the following information.
 *************
		*****************************************************************
		1)The first line: vertexNum	picking_ratio	degree.
		*****************************************************************
		2)The second part is the task assignment graph, GT.
		*****************************************************************
		3)The third part is the Ground truth ranking: For example, 0 > 1 > 2 > 3 > 4 > 5 > 6 > 7 > 8 > 9 > 10.
		*****************************************************************
		4)The forth part is the voting matrix.
		*****************************************************************
		5)The fifth part is preference graph GP.
		*****************************************************************
		6)The sixth part is transitivity closure graph GTC.
		*****************************************************************
		7)This part is the result your algorithm (i.e., truth discovery) needs. 
		It arranges as follows: one line for one worker, consisting of 4-tuple data, that is,
        <1-st is task-id, 2-nd is object i, 3-rd is object j, 4-th is the preference result>.
        Please note that all the indices are zero-based. 
*/

#define _GT_VM_GP_Generator_CPP_ /*comment this line, or*/
#ifndef _GT_VM_GP_Generator_CPP_ /*change ifndef to ifdef, to make this main function work.*/

#define _CRT_SECURE_NO_WARNINGS
//#define __RELEASE_MODE_CORWDSOURCING_PROJECT_

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
#include "truth_discovery.hpp"

using namespace std;


int main(int argc, char * argv[]){

	string s_vertexNum = argv[1]; // argv[3];
	string s_ratio = argv[2]; // argv[4];
	string s_Gauss_mean = argv[3];
	// Gaussian variance to control the worker's quality.
	string s_Gauss_stddev = argv[4];
	string s_distribuion_type = argv[5]; // argv[9]; == 1(Gauss); == 0(Uniform);
	string s_workerNum = argv[6];
	string file_path = argv[7];

	

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
		cout << "GT cannot be generated. Try again!\n";
		return -1;
	}

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string(now->tm_hour)
		+ "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec);

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

	// generate voting matrix, VM in short.
	// ground truth is : 0 > 1 > 2 > ... > n-1.
	int * GTr = new int[vertexNum];
	std::map<int, int> invert_index;

	of1 << "\nGround truth is : \n";
	for (int i = 0; i < vertexNum; i++){	
		GTr[i] = i;
		invert_index[GTr[i]] = i;
		of1 << i << " > ";
	}
	of1 << "\n";

	int edge_num = get_edge_num(g_t, vertexNum);
	// maximum number of pairwise comparison per worker can do.
	double tasks_ratio = 0.01;
	int max_per_worker = floor((double)edge_num * tasks_ratio);
	if (max_per_worker < 3) max_per_worker = 3;
	int workerNum = stoi(s_workerNum);
	int distribuion_type = stoi(s_distribuion_type);
	double mean = stod(s_Gauss_mean);
	double stddev = stod(s_Gauss_stddev);
	std::srand(unsigned(std::time(0)));
	vector<vector<vector<int>>> matrice= generate_sim_votes(vertexNum, workerNum,
		g_t, invert_index, distribuion_type, mean, stddev, max_per_worker);

	int smothingType = 3;/*LS, i.e., Laplace smoothing.*/
	std::vector<std::vector<double>> gp = 
		generate_aggregated_matrix(matrice, smothingType);
	std::vector<std::vector<double>> gp_tc = generateTransitionMatrix_Eigen(gp);
	of1 << "\n-------------------------\n";
	of1 << "\nVoting matrix is: \n";
	for (int ww = 0; ww < workerNum; ++ww){
		of1 << "Worker " << ww << " has the results:\n";
		fprint_2D_vector<int>(vertexNum, matrice[ww], of1);
	}

	of1 << "\n-------------------------\n";
	of1 << "\nGP is: \n";
	fprint_2D_vector<double>(vertexNum, gp, of1);

	of1 << "\n-------------------------\n";
	of1 << "\nGTC is: \n";
	fprint_2D_vector<double>(vertexNum, gp_tc, of1);

	vector<vector <td_wieghts::client_vote>> v_client_votes =
		vector<vector<td_wieghts::client_vote>>(workerNum, vector<td_wieghts::client_vote>(0));
	generate_sim_votes_taskID_who(vertexNum, workerNum, g_t, matrice, v_client_votes);


	// store the tasks;
	vector <pair<int, int>> v_tasks;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < i; ++j){
			if (g_t[i][j] == true){
				v_tasks.push_back(make_pair(i, j));
			}
		}
	}

	// ??????????????
	// each work will vote for some tasks;
	of1 << "\n-------------------------\n";
	of1 << "The result for truth discovery method.\n";
	for (int ww = 0; ww < workerNum; ++ww){
		for (auto & task : v_client_votes[ww]){
			of1 << task.task_id << "\t" << v_tasks[task.task_id].first <<
				"\t" << v_tasks[task.task_id].second <<
				"\t" << task.who << "\t";
		}
		of1 << "\n";
	}

	of1.close();
	// release memory;
	delete_2D_Array<bool>(vertexNum, g_t);
	delete[]  GTr;
	std::vector<std::vector<double>>().swap(gp);
	std::vector<std::vector<double>>().swap(gp_tc);
	vector<vector<vector<int>>>().swap(matrice);
	vector<vector <td_wieghts::client_vote>>().swap(v_client_votes);
	return 0;
}
#endif