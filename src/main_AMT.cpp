/**
* @file: main_AMT.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#define _MAIN_AMT__FUNCTION_CPP_ /*comment this line, or*/
#ifndef _MAIN_AMT__FUNCTION_CPP_ /*change ifndef to ifdef, to make this main function work.*/
#define _CRT_SECURE_NO_WARNINGS

#include "quickSort.hpp"
#include "MyUtility.hpp"
#include "ta.hpp"
#include "myGraph.hpp"
#include "readCSVFiles.hpp"
#include "greedyDFS.hpp"
#include "simulation.hpp"
#include "printParams.hpp"
#include <opencv2/highgui/highgui.hpp>
#include "ktaub.hpp"
#include <ctime>
#include "heurisHP.hpp" /*class NormalFormMatrix4HP*/
#include "graph_ta.hpp"
#include "makeDirectory.hpp"
#include "TSPalgorithm.hpp"
#include "matrix_generator.hpp"
#include "transition.hpp"
#include <chrono> /*time*/
#include "Experiment.hpp"
#include "Baseline2.hpp"
#include "repeatchoice.hpp"
#include <fstream>
#include <iomanip> /*IO Manipulators*/
#include "task_assign.hpp"
#include "CrowdBT.hpp"
#include <cstdio> /*printf*/
#include "truth_discovery.hpp"
#include "ta.hpp"


// #define __RELEASE_MODE_CORWDSOURCING_PROJECT_
using namespace cv;
using namespace std;
extern std::default_random_engine generator_main;
extern int  totalTime;
extern double beta_a, beta_b, _ACCURCY;
//double beta_a = 10;
//double beta_b = 5;
//double _ACCURCY = 1.0E-10;
//int  totalTime = 50000;


using namespace std;


int main(int argc, char * argv[]){
	// seed for rand(), shuffle_rand(), etc.
	std::srand(unsigned(std::time(0)));
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	generator_main = std::default_random_engine(seed);

	const double coolrate = 0.95;
	const int IterationNum = 1000;
	// temperature;
	const int temperature = 50000;
	const int ExperimentTimes = 10;

	string nam = "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/AMT/SATD/workers_quality.txt";
	std::ofstream  of(nam.c_str(), std::ofstream::out | std::ofstream::app);

	FILE *fin;
	vector<string> gtfilelist;
	string baseAddress = "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/AMT/G_P_For_Each_Worker/together/";
	string outputAddress = "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/AMT/SATD/together/";
	GetFileList(baseAddress, ".txt", &gtfilelist);


	// each GT file;
	for (int fidx = 0, filesSize = gtfilelist.size(); fidx < filesSize; fidx++){
		cout << "Processing " << gtfilelist[fidx] << ":\n";
		cout << " --- Reading voting matrix for each worker: " << gtfilelist[fidx] << endl;
		// read GT file
		string filename = baseAddress + gtfilelist[fidx];
		int vertexNum = 0, workerNum = 0;
		fin = fopen(filename.c_str(), "r");
		// "%lf" for double format.
		fscanf(fin, "%d\t%d", &vertexNum, &workerNum);
		cout << " vertexNum = " << vertexNum << ", workerNum = " << workerNum << endl;

		// new GT
		bool ** g_t = build_2D_array<bool>(vertexNum, false);
		double ** g_tc = build_2D_array<double>(vertexNum, 0.0);

		vector<vector<vector<int>>> matrice;
		vector<vector<int>> temp_voting(vertexNum, vector<int>(vertexNum, 0));

		char *buff = (char*)malloc(sizeof(char) * 4096);
		int client_n = 0;
		while (1){
			// reads a line from the specified stream and stores it into the string;
			// It stops when either (n-1) characters are read, 
			// the newline character is read, or the end-of-file is reached, 
			// whichever comes first.
			fgets(buff, 4096, fin);
			// Preference graph G_P for C++:
			if (buff[0] == 'P' && buff[1] == 'r' && buff[2] == 'e' && buff[3] == 'f'){
				client_n++;
				for (int i = 0; i < vertexNum; ++i){
					fgets(buff, 4096, fin);
					for (int k = 0, j = 0; j < vertexNum; ++j, k += 3){
						//cout << "buff[" << k << "] = " << buff[k] << ", ";
						temp_voting[i][j] = buff[k] - '0';
					}
				}
				//cout << "Worker " << client_n << ":\n";
				//print_2D_vector<int>(vertexNum, temp_voting);
				matrice.push_back(temp_voting);
			}
			if (client_n == workerNum)
				break;
		}

		/*for (int i = 0; i < workerNum; ++i){
			cout << "worker " << i << " votes:\n";
			print_2D_vector<int>(vertexNum, matrice[i]);
		}*/

		td_wieghts::truth_discovery * td = new td_wieghts::truth_discovery();
		int max_task_num = 1000;
		td->init(1, vertexNum, workerNum, max_task_num);
		int smothingType = 1;
		int TruthIterationTimes = 0;
		std::vector<std::vector<double>>
			gp_SATD_wwm = generate_aggregated_matrix(matrice, smothingType);
		int edge_num = 0;
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < i; ++j){
				if (gp_SATD_wwm[i][j] + gp_SATD_wwm[j][i] > 0){
					g_t[i][j] = true;
					g_t[j][i] = true;
					edge_num++;
				}
				else{
					g_t[i][j] = false;
					g_t[j][i] = false;
				}
			}
		}
		vector<pair<int, int>> e_deges;
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < i; ++j){
				if (g_t[i][j] || g_t[j][i])
					e_deges.push_back(make_pair(i, j));
			}
		}
		int temp = 0;

		/*for (auto i : e_deges){
			
			cout << "edge " << temp << ": " << i.first << ", " <<
				i.second << endl;
				temp++;
		}*/

		//print_2D<bool>(vertexNum, g_t);
		//print_2D_vector<double>(vertexNum, gp_SATD_wwm);

		if (!gp_SATD_wwm.empty())
			gp_SATD_wwm.clear();

		vector<double> worker_weights(workerNum, 0.0);
		vector<vector <td_wieghts::client_vote>> v_client_votes =
			generate_aggregated_matrix_SATD(g_t, matrice, gp_SATD_wwm,
			worker_weights, td, max_task_num, totalTime, smothingType);

		
		double prob_threshold = 0.01;
		double stddev_scalar = 1.0;
		double mean = 0;
		Smoothing_UpdateProb_with_Workers_Quality(
			worker_weights, g_t, vertexNum, prob_threshold,
			mean, stddev_scalar, v_client_votes,
			gp_SATD_wwm);
		//print_2D_vector<double>(vertexNum, gp_SATD_wwm);
		std::vector<std::vector<double>>
			gp_tc_SATD_wwm = generateTransitionMatrix_Eigen(gp_SATD_wwm);

		//print_2D_vector<double>(vertexNum, gp_tc_SATD_wwm);

		string tem_name = outputAddress + gtfilelist[fidx];
		std::ofstream  of2(tem_name.c_str(), std::ofstream::out);
		
		int E = e_deges.size();
		vector<double > v_weights(vertexNum*vertexNum, 0.0);
		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				v_weights[i*vertexNum + j] = gp_tc_SATD_wwm[i][j];

		vector<int> vertexSet;
		for (int i = 0; i < vertexNum; i++){
			vertexSet.push_back(i);
		}

		int * TAr = new int[vertexNum];
		int * TAr_inv = new int[vertexNum];

		TA ta(vertexNum, E, 10 /*top_k*/, false, false);
		ta.initiaWeights(v_weights);
		ta.bruteForce(vertexSet, of2);
		TAr = ta.getmaxPath();

		for (int i = 0; i < vertexNum; i++){
			TAr_inv[i] = TAr[vertexNum - i - 1];
			std::cout << TAr[i] << ", ";
			of2 << TAr[i] << ", ";
		}

		if (!of2.is_open())
			cout << "Not Open file " << tem_name << ".\n";


		//Run SATD_TT From TC
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				g_tc[i][j] = gp_tc_SATD_wwm[i][j];
			}
		}

		int d = 2 * edge_num / vertexNum;
		int permuStartingIdx = 1;
		int distribuion_type = 1;
		int SA_flag = 2;
		double sum_ = 0;
		for (int rep = 0; rep < ExperimentTimes; ++rep){
			cout << gtfilelist[fidx] << " repeat " << rep << "/" << ExperimentTimes
				<< " trails.\n";
			double tau = 0, sec = .0;
			Run_SA_SATD_From_TC(tau, sec,
				TAr, // groud truth ranking;
				matrice, g_tc,
				vertexNum, d, coolrate, IterationNum, temperature,
				permuStartingIdx, distribuion_type,
				workerNum, SA_flag, "SATD_WWM", of2);
			sum_ += tau;

			of2 << gtfilelist[fidx] <<  " Trial " << rep << ": tau = " << tau << endl;

		}/*end of each time experiment*/

		of2 << gtfilelist[fidx] << " average tau = " << sum_/ ExperimentTimes << endl;

		std::size_t found2 = gtfilelist[fidx].find("_Ea");

		of << string(gtfilelist[fidx], 0, found2) << "Execution Iteration times is " 
			<< TruthIterationTimes << " / " << totalTime << ", Workers' quality are:\n";
		for (auto i : worker_weights)
			of << i << "\, ";
		of << "\n\n";

		fclose(fin);
		of2.close();
		// release memory;
		delete_2D_Array<bool>(vertexNum, g_t);
		delete_2D_Array<double>(vertexNum, g_tc);
		td->free_all();
		delete td;
	}/*end of each GT file*/
	return 0;
}
#endif