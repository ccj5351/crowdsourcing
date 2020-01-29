/**
* @file: Experiment.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_EXPERIMENTS_H_
#define __HEADER__CROWD_SOURCING_EXPERIMENTS_H_

#include "MyUtility.hpp"
#include "ta.hpp"
#include "myGraph.hpp"
#include "readCSVFiles.hpp"
#include "greedyDFS.hpp"
#include "simulation.hpp"
#include "printParams.hpp"
#include <opencv2/highgui/highgui.hpp>
#include "factorial2int.hpp"
#include "ktaub.hpp"
#include <ctime>
#include "heurisHP.hpp" /*class NormalFormMatrix4HP*/
#include "graph_ta.hpp"
#include "makeDirectory.hpp"
#include "TSPalgorithm.hpp"
#include "matrix_generator.hpp"
#include "transition.hpp"
#include "repeatchoice.hpp"
#include <chrono> /*time*/
#include <algorithm>

#define _TURE_OR_FALSE_ 1
#define THRESHOLD_TIME_IN_SECONDS 180 // parameters to control the time spent on task assignment.
#define TASK_ASSIGNMENT_MAX_TIMES 10  // parameters to control the time spent on task assignment.
using namespace std; 

void TA_SA_Compare(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	const int & d,// degree
	const int & parallelIdx,
	const double &coolrate, // e.g. = (float) 0.95;
	const int & IterationNum, // e.g., = 100000;
	const int & temperature, // temperature, e.g., = 100000;
	const int & ExperimentTimes, // how many random trials will be  
	// executed, for an averaged result.
	const int & permuStartingIdx, // e.g., = 1. means that for any permutation,
	// the first (i.e., 0) vertex will not be moved.
	const int & distribuion_type,
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & max_per_worker, // maximum number of pairwise comparison per worker can do.
	/*Smothing_Type : */
	const int & smothingType,
	const int & SA_flag,
	const int & workerNum,
	std::ofstream & of1
	);

void TA_GroundTruth_Compare(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	const int & d,// degree
	const int & parallelIdx,
	const int & ExperimentTimes, // how many random trials will be  
	// executed, for an averaged result.
	const int & distribuion_type,
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1 
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & max_per_worker, // maximum number of pairwise comparison per worker can do.
	/*Smothing_Type : */
	const int & smothingType,
	std::ofstream & of1
	);

// try to do task assignment to generate G_T graph;
// if the degree d is not appropriate,
// we will try some other degree within some time threshold;
// otherwise, false will be returned.
bool task_assignment_graph(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	int & d// degree 
	// ofstream & of1
	);

bool get_G_T(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	const double & thre, // time threshold, since we do not
	// want the task assignment process falls into 
	// an endless loop.
	const int & d// degree 
	);


int SA_SATD_CBT_RC_GM_Comparison(
	double & SATD_tau, // SA with truth discovery tau distance;
	double & CBT_tau, // CBT with truth discovery tau distance;
	double & RC_tau, // RC tau distance;
	double & GM_tau, // GM tau distance;
	double & SATD_sec, // SA with truth discovery running time;
	double & CBT_sec, // CBT running time;
	double & RC_sec, // RC running time;
	double & GM_sec, // GM running time;
	double & SATD_WWM_tau, // SATD_WWM tau distance;
	double & SATD_WWM_sec, // SATD_WWM total running time;
	double & SATD_WWM_TD_sec, // SATD_WWM running time, for truth discovery.
	double & SATD_WWM_Smoothing_sec, // SATD_WWM running time, for smoothing.
	double & SATD_WWM_TC_sec, // SATD_WWM running time, for Transitive closure.
	double & SATD_WWM_HP_sec, // SATD_WWM running time, for finding best HP.
	const int & vertexNum, // number of vertex ;
	int & d,// degree
	bool ** g_t, // task assignment graph.
	const double & ratio,
	const double & coolrate, // e.g. = (float) 0.95;
	const int & IterationNum, // e.g., = 100000;
	const int & temperature, // temperature, e.g., = 100000;
	const int & ExperimentTimes, // how many random trials will be  
	// executed, for an averaged result.
	const int & CBT_ExperimentTimes, // how many random trials will be  
	// executed, for Crowd_BT, since it consumes too much time.
	const int & permuStartingIdx, // e.g., = 1. means that for any permutation,
	// the first (i.e., 0) vertex will not be moved.
	const int & distribuion_type,
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1 
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & smothingType, /*Smothing_Type : */
	const int & SA_flag,
	const int & workerNum,
	const double & tasks_ratio,
	//const bool & IsTrustTransitivityFlag, // == true, using the TrustTransitivity, equation (7).
	// weight of direct edge probability.
	// const double & alpha, // alpha* prob[i][j] + (1 - alpha)*a[i][j];
	/*a threshold of the small_prob to be generated.*/
	const double & prob_threshold, 
	const bool & Is_Beta_Simu_data_type,
	vector <int> & TruthIterationTimes,
	std::ofstream & of1
	);


void Run_SA_SATD_From_TC(
	double & tau,
	double & sec,
	// the same order as SA, i --> j, means to prefer j than i.
	const int * GTr, // groud truth ranking;
	const std::vector<std::vector<std::vector<int>>> & v_M, // each voter's voting result;
	double ** p_Gp_tc, // transitive closure matrix G_TC;
	const int & vertexNum, // number of vertex ;
	const int & d,// degree
	const double &coolrate, // e.g. = (float) 0.95;
	const int & IterationNum, // e.g., = 100000;
	const int & temperature, // temperature, e.g., = 100000;
	const int & permuStartingIdx, // e.g., = 1. means that for any permutation,
	// the first (i.e., 0) vertex will not be moved.
	const int & distribuion_type,
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1 
	const int & SA_flag,
	const int & workerNum,
	const string & TaskName, // e.g., = SA, or  = SATD;
	std::ofstream & of1
	);

void Run_RC_From_VotingMatrix(
	double & tau,
	double & sec,
	// pay attention here:
	// the order is i --> j, means to prefer i than j, is different from SA;.
	const int * GTr, // groud truth ranking;
	const vector<vector<vector<int> >> & matrice,
	const int & vertexNum, // number of vertex ;
	const int & d,// degree
	const int & ExperimentTimes, // how many random trials will be  
	const int & distribuion_type,
	const int & w, // number of workers;
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1 
	std::ofstream & of1
	);

void Run_GM_From_VotingMatrix(
	double & tau,
	double & sec,
	// pay attention here:
	// the order is i --> j, means to prefer i than j, is different from SA;.
	const int * GTr, // groud truth ranking;
	const vector<vector<vector<int> >> & matrice,
	const int & vertexNum, // number of vertex ;
	const int & d,// degree
	const int & distribuion_type,
	const int & w, // number of workers;
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1 
	std::ofstream & of1
	);


template <typename T>
T **  build_2D_array(const int & vertexNum, const T & initial){
	T ** p_Gp_tc = new T * [vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		p_Gp_tc[i] = new T [vertexNum];
	}
	for (int i = 0; i < vertexNum; ++i)
		for (int j = 0; j < vertexNum; ++j)
			p_Gp_tc[i][j] = initial;
	return p_Gp_tc;
}


template <typename T>
T  ***  build_3D_array(const int & workerNum, const int & vertexNum, const T & initial){
	T *** p_matrice = new T **[workerNum];
	for (int i = 0; i < workerNum; i++){
		p_matrice[i] = new T *[vertexNum];
		for (int j = 0; j < vertexNum; j++)
			p_matrice[i][j] = new T[vertexNum];
	}

	for (int i = 0; i < vertexNum; ++i)
		for (int j = 0; j < vertexNum; ++j)
			for (int k = 0; k < vertexNum; ++k)
				p_matrice[i][j][k] = initial;
	return p_matrice;
}

template <typename T>
void delete_3D_Array(const int & w, const int & n, T *** p_matrice){
	for (int i = 0; i < w; i++){
		for (int j = 0; j < n; i++){
			delete[] p_matrice[i][j];
		}
	}
	delete[] p_matrice;
}


template <typename T>
void delete_2D_Array(const int& n, T **p_2D) {
	for (int i = 0; i < n; i++)
		delete[] p_2D[i];
	delete[] p_2D;
}

template<typename T>
void print_2D(const int& n, T ** pp){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			cout << pp[i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

template<typename T>
void fprint_2D(const int& n, T ** pp, ofstream & of){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			of << pp[i][j] << ", ";
		of << endl;
	}
	of << endl;
}

template<typename T>
void print_2D_vector(const int& n, const vector<vector<T>> & pp){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			cout << pp[i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}


inline void print_2D_Eigen(const int& n, const Eigen::MatrixXd & pp){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			cout << pp(i,j) << ", ";
		cout << endl;
	}
	cout << endl;
}

template<typename T>
void fprint_2D_vector(const int& n, const vector<vector<T>> & pp,
	ofstream & of){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			of << pp[i][j] << ", ";
		of << endl;
	}
	of << endl;
}

template<typename T>
void print_1D(const int& n, T * p){
	for (int i = 0; i < n; ++i)
		cout << p[i] << ", ";
	cout << endl;
}

void set_GT_False(const int& n, bool **GT);


std::vector<std::vector<double>> getGTCfromGPInFile(const string path,
	const int & vertexNum);

bool getGTCfromGPInFile(const string path, double ** p_gtc,
	const int & vertexNum);

void crowdsourcing_Compare_TA_SA(
	const string & baseAddressGT,
	const string & baseAddressGP,
	const string & baseAddressTC,
	const int & vertexNum,
	const int & d,
	const int & IterationNum,
	const int & permuStartingIdx,
	const int & temperature,
	const double & coolrate,
	const int & ExperimentTimes,
	const int & SA_flag,
	const int & w // number of workers;
	);

inline void printParameters(
	const string & title,
	const int & vertexNum,
	const int & d,
	const int & temperature,
	const double & coolrate,
	const int &IterationNum,
	const int & ExperimentTimes,
	const int & distribuion_type,
	const int & workerNum,
	const int & permuStartingIdx,
	const bool & isDisplay,
	const string & cur_time
	){

	std::cout << " //*************************************" << endl;
	std::cout << " //******" << title << "*********" << endl;
	std::cout << " //*************************************" << endl;
	std::cout
		<< "Parameters include:\n"
		<< "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl
		<< "Initial Temperature T = " << temperature << endl
		<< "Cooling rate = " << coolrate << endl
		<< "Iteration Number = " << IterationNum << endl
		<< "Experiment Times for average result = " << ExperimentTimes << endl;
	if (distribuion_type == 0)
		cout << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		cout << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	cout << "workerNum = " << workerNum << endl;
	cout << "permuStartingIdx = " << permuStartingIdx << endl;
	//cout << "Parallel vertex index if needed = " << parallelIdx << endl;
	if (isDisplay)
		cout << "isDisplay = True" << endl;
	else
		cout << "isDisplay = False" << endl;
	cout << "Current time is " << cur_time << endl;
}

inline void saveParameters(
	const string & title,
	const int & vertexNum,
	const int & d,
	const int & temperature,
	const int & coolrate,
	const int &IterationNum,
	const int & ExperimentTimes,
	const int & distribuion_type,
	const int & workerNum,
	const int & permuStartingIdx,
	const bool & isDisplay,
	const string & s_tasks_ratio,
	const string & cur_time,
	ofstream & of
	){
	of << " //*************************************" << endl;
	of << " //******" << title << "*********" << endl;
	of << " //*************************************" << endl;
	of
		<< "Parameters include:\n"
		<< "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl
		<< "Initial Temperature T = " << temperature << endl
		<< "Cooling rate = " << coolrate << endl
		<< "Iteration Number = " << IterationNum << endl
		<< "Experiment Times for average result = " << ExperimentTimes << endl;
	if (distribuion_type == 0)
		of << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		of << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	of << "workerNum = " << workerNum << endl;
	of << "permuStartingIdx = " << permuStartingIdx << endl;
	//of << "Parallel vertex index if needed = " << parallelIdx << endl;
	if (isDisplay)
		of << "isDisplay = True" << endl;
	else
		of << "isDisplay = False" << endl;
	of << "Tasks ratio is " << s_tasks_ratio << endl;
	of  << "Current time is " << cur_time << endl;
}

inline void saveParameters(
	const string & title,
	const int & vertexNum,
	const int & d,
	const int & temperature,
	const int & coolrate,
	const int &IterationNum,
	const int & ExperimentTimes,
	const int & distribuion_type,
	const int & workerNum,
	const int & permuStartingIdx,
	const bool & isDisplay,
	const string & cur_time,
	ofstream & of
	){
	of << " //*************************************" << endl;
	of << " //******" << title << "*********" << endl;
	of << " //*************************************" << endl;
	of
		<< "Parameters include:\n"
		<< "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl
		<< "Initial Temperature T = " << temperature << endl
		<< "Cooling rate = " << coolrate << endl
		<< "Iteration Number = " << IterationNum << endl
		<< "Experiment Times for average result = " << ExperimentTimes << endl;
	if (distribuion_type == 0)
		of << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		of << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	of << "workerNum = " << workerNum << endl;
	of << "permuStartingIdx = " << permuStartingIdx << endl;
	//of << "Parallel vertex index if needed = " << parallelIdx << endl;
	if (isDisplay)
		of << "isDisplay = True" << endl;
	else
		of << "isDisplay = False" << endl;
	of << "Current time is " << cur_time << endl;
}

inline void saveParameters(
	const string & title,
	const int & vertexNum,
	const int & d,
	const int & distribuion_type,
	const int & workerNum,
	const string & cur_time,
	ofstream & of
	){
	of << " //*************************************" << endl;
	of << " //******" << title << "*********" << endl;
	of << " //*************************************" << endl;
	of
		<< "Parameters include:\n"
		<< "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl;
	if (distribuion_type == 0)
		of << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		of << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	of << "workerNum = " << workerNum << endl;
	of << "Current time is " << cur_time << endl;
}

inline std::string generateFileName(
	const string & path, // "path"
	const int & vertexNum,
	const int & d,
	const int & IterationNum,
	const int & distribuion_type,
	const int & SA_flag,
	const int & workerNum,
	const double & ratio,
	const string & s_task_name,
	const string & cur_time,
	const string & file_extension // e.g., = ".txt"
	){
	string temp;
	string s_flag;
	if (distribuion_type == 1)
		temp = "_Gaus";
	else if(distribuion_type == 0)
		temp = "_Unif";
	else
		temp = "_Other";

	if (SA_flag == 0)
		s_flag = "_RTour_";
	else if (SA_flag == 1)
		s_flag = "_NNTour_";
	else if (SA_flag == 2)
		s_flag = "_RankVTour_";
	else if (SA_flag == 3)
		s_flag = "_RankTour_";
	else if (SA_flag == 4)
		s_flag = "_CondorFusTour_";
	else
		s_flag = "_NNTour_";

	string temp1 = to_string(log10(IterationNum));
	string temp2 = to_string(ratio);
	string tem_name = path + "V_" + to_string(vertexNum) + "_d_" + to_string(d) + "_W_" + to_string(workerNum)
	+ "_N_" + string(temp1.begin(), temp1.begin()+ 4 ) + temp + s_flag + "r_" +
	string(temp2.begin(), temp2.begin() + 4) + "_" + s_task_name  // = "_TA_SA_Result.txt"
		+ "_" + cur_time + file_extension;
	return tem_name;
}

bool Run_Make_GT(
	const int & vertexNum,
	int & d,
	bool ** g_t,
	ofstream & of
	);

inline void printInitalTourFinalTour(
	const double & iniScore,
	const double & score,
	const bool & iniIsHP,
	const bool & ResultIsHP
	){

	std::cout << "\ninitialScore/ResultScore = " << iniScore << "/" << score
		<< ", iniIsHP/ResultIsHP = ";
	if (iniIsHP)
		cout << "True/";
	else cout << "False/";

	if (ResultIsHP)
		cout << "True\n";
	else cout << "False\n";
}

// get the number of edges (i.e., the number of pairwise comparison tasks).
// GT: task assignment graph.
inline int get_edge_num(bool ** GT, const int & vertexNum){
	int edge_num = 0;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = i + 1; j < vertexNum; ++j){
			if (GT[i][j] == true)
				edge_num++;
		}
	}
	return edge_num;
}

void run_CrowdBT_for_SA_Comparison(
	double & tau,
	double & sec,
	const int & vertexNum,
	// pay attention here:
	// the order is i --> j, means to prefer i than j, is different from SA;.
	const int * GTr, // groud truth ranking; 
	const int & distribuion_type,
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & w, // number of workers;
	const int & num_repeat,
	const bool & Is_Beta_Simu_data_type // == true, using Beta distribution to simulate the voting result.
	);

#endif