/**
* @file: main_SATD.cpp
* @brief: The main function for the experiments of our paper.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#define _MAIN_main_SA_SATD_CrowdBT_CPP_ /*comment this line, or*/
#ifdef  _MAIN_main_SA_SATD_CrowdBT_CPP_ /*change ifndef to ifdef, to make this main function work.*/
#define _CRT_SECURE_NO_WARNINGS
#define __RELEASE_MODE_CORWDSOURCING_PROJECT_

#include "quickSort.hpp"
#include "MyUtility.hpp"
#include "ta.hpp"
#include "myGraph.hpp"
#include "readCSVFiles.hpp"
#include "greedyDFS.hpp"
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


// #define __RELEASE_MODE_CORWDSOURCING_PROJECT_
using namespace cv;
using namespace std;
std::default_random_engine generator_main;
vector<double> v_td_secs; // truth discovery;
vector<double> v_tc_secs; // transitive
vector<int> v_smoothing_counts; 
vector<double> v_smoothing_time_debug(1000, 0);
int temp_main;

double beta_a = 10;
double beta_b = 5;
double _ACCURCY = 1.0E-10;
int  totalTime = 6000;
double ratio_main;


/*a threshold of the small_prob to be generated.*/
// double prob_threshold;
#if 1


int main(int argc, char * argv[]){

	string program = argv[0];
	string baseAddress = argv[1];
	string sk = argv[2]; // finding top-k paths using TA;  //argv[2];
	string s_vertexNum = argv[3]; // argv[3];
	string s_d = argv[4]; // argv[4];
	string s_T = argv[5]; // argv[5];
	string s_coolRate = argv[6]; // argv[6];
	string s_IterationNum = argv[7]; // argv[7];
	string s_ExperimentTimes = argv[8]; // argv[8];
	string s_distribuion_type = argv[9]; // argv[9];
	string s_permuStartingIdx = argv[10]; // argv[10];
	// number of workers for each pairwise comparison, 
	// that is, how many workers will work on the pairwise
	// comparison (i.e., one edge).
	string s_workerNum = argv[11];
	string s_ratio = argv[12];
	string s_SA_flag = argv[13];
	int arg_idx = 14;
	string s_isDisplay = argv[arg_idx];
	bool isDisplay = s_isDisplay == "1" ? true : false;
	bool is_Debug = false;
	string s_tasks_ratio = argv[arg_idx + 1];
	// different levels of error-rate for the worker's quality.
	string s_Worker_quality = argv[arg_idx + 2];
	// Gaussian variance to control the worker's quality.
	// string s_Gauss_stddev = argv[arg_idx + 3];
	string resutsDir = argv[arg_idx + 3];
	string s_smothingType = argv[arg_idx + 4];
	
	// weight of direct edge probability.
	// alpha* prob[i][j] + (1 - alpha)* a[i][j];
	// where prob[i][j] means direct edge probability between node i and j;
	// and a[i][j] means the inferred (i.e., the transitivity closure) 
	// probability between node i and j;
	string s_prob_threshold = argv[arg_idx + 5];
	/*a threshold of the small_prob to be generated.*/
	double prob_threshold = stod(s_prob_threshold);
	string s_Is_Beta_Simu_data_type = argv[arg_idx + 6];
	bool Is_Beta_Simu_data_type = (stoi(s_Is_Beta_Simu_data_type) == 1) ?
		true: false;

	// e.g., s_task_name ="_GreedyDFS", or "_TA", etc.;
	string s_task_name = argv[arg_idx + 7];
	
	// seed for rand(), shuffle_rand(), etc.
	std::srand(unsigned(std::time(0)));

	int vertexNum = stoi(s_vertexNum);
	int K = stoi(sk);
	//int parallelIdx = stoi(s_parallel);
	int permuStartingIdx = stoi(s_permuStartingIdx);
	int SA_flag = stoi(s_SA_flag);

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	/*string cur_time = to_string((now->tm_mon + 1)) + "-" 
		+ to_string(now->tm_mday) + "-" + to_string(now->tm_hour) 
		+ "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec);*/
	string cur_time = to_string(now->tm_hour)
		+ "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec);

	const double coolrate = stod(s_coolRate);// e.g. = (float) 0.95;
	const int IterationNum = stoi(s_IterationNum);// e.g., = 100000;
	// temperature;
	const int temperature = stoi(s_T);// e.g.,100000;
	const int ExperimentTimes = stoi(s_ExperimentTimes);// e.g.,20;
	// #define UNIFORM_DISTRIBUTION 0
	// #define NORMAL_DISTRIBUTION 1
	const int distribuion_type = stoi(s_distribuion_type);
	int workerNum = stoi(s_workerNum);
	
	v_tc_secs = vector<double>(ExperimentTimes, 0);
	v_td_secs = vector<double>(ExperimentTimes, 0); // truth discovery;
	v_tc_secs = vector<double>(ExperimentTimes, 0); // transitive
	v_smoothing_counts = vector<int>(ExperimentTimes, 0);
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	generator_main = std::default_random_engine(seed);

	// new change, paper for ICDE 2017;
	// run SA, RC, GM. Crowd_BT together, from the same simulation data;
	double ratio = stod(s_ratio);
	double tasks_ratio = stod(s_tasks_ratio);
	int d = 0;

	string title = "Comparing SATD, SATD_WS, RC, GM, and Crowd_BT.";
	printParameters(title, vertexNum, d, temperature,
		coolrate, IterationNum, ExperimentTimes, distribuion_type,
		workerNum, permuStartingIdx, isDisplay, cur_time);

	// saveing result;
	string tem_name;
	if (distribuion_type == 0)
		tem_name = resutsDir + "v" + s_vertexNum + "-w" + s_workerNum 
		+ "-rt-" + s_tasks_ratio + "-Uniform-WQ-" + s_Worker_quality 
		+ "-Smooth-" + s_smothingType + "-Time-" + cur_time + ".txt";
	if (distribuion_type == 1)
		tem_name = resutsDir + "v" + s_vertexNum + "-w" + s_workerNum +
		"-rt-" + s_tasks_ratio + "-Gauss-WQ-" + s_Worker_quality 
		+ "-Smooth-" + s_smothingType + "-Time-" + cur_time + ".txt";
	std::ofstream of1(tem_name.c_str(), std::ofstream::out);

	if (!of1.is_open())
		cout << "Not Open file " << tem_name << ".\n";


	FILE *fin;
	vector<string> gtfilelist;
	GetFileList(baseAddress, ".txt", &gtfilelist);
	int CBT_ExperimentTimes = 1;
	
	int smothingType = stoi(s_smothingType);/*Smothing_Type : */
	// Gaussian variance to control the worker's quality.
	int Worker_quality = stoi(s_Worker_quality);
	// save the task assignment graph.
	of1 << "#-----------------------------------\n";
	of1 << "# vertexNum = " << vertexNum << "\tworkerNum = "
		<< workerNum << "\tRatio = "
		<< ratio << "\tDegree = " << d << "\tDistribution = " << distribuion_type
		<< "\n#Worker Quality Type = " << Worker_quality
		<< "\tExperimentTimes = " << ExperimentTimes
		<< "\tIterationNum = " << IterationNum
		<< "\tBeta_a = " << beta_a
		<< "\tBeta_b = " << beta_b
		<< "\tSmoothing Type = " << smothingType
		<< "\tTask Ratio = " << tasks_ratio 
		<< "\n#Is_Beta_Simu_data_type = " << Is_Beta_Simu_data_type
		<< "\tSmoothingProbThreshold = " << prob_threshold
		<< "\tCBT_ExperimentTimes = " << CBT_ExperimentTimes
		<< "\tSA_flag = " << SA_flag << "\n";

	of1 << "#ratio\t" 
		<< "SATD_tau\t"  << "SATD_WWM_tau\t"   << "RC_tau\t"  << "GM_tau\t"  << "CBT_tau\t"
		<< "SATD_time\t" << "SATD_WWM_time\t"  << "RC_time\t" << "GM_time\t" << "CBT_time\t"
		<< "VertexNum\t" << "Seperate_TD_time\t" << "Seperate_Smoothing_time\t"
		<< "Seperate_TC_time\t" << "Seperate_HP_time\n";


	   // SATD_WWM running time, for truth discovery.
	   // SATD_WWM running time, for smoothing.
       // SATD_WWM running time, for Transitive closure.
	   // SATD_WWM running time, for finding best HP.

	for(int fidx = 0, filesSize = gtfilelist.size();
		fidx < filesSize; fidx++){ // each GT file;
		cout << "\nProcessing GT : " << gtfilelist[fidx] << endl;
		// read GT file
		string filename = baseAddress + gtfilelist[fidx];
		fin = fopen(filename.c_str(), "r");
		// "%lf" for double format.
		fscanf(fin, "%d\t%lf\t%d", &vertexNum, &ratio, &d);
		
		if (d == vertexNum - 1)
			ratio = 1.0;
		else
			ratio = double(d) / double(vertexNum);

		ratio_main = ratio;
		// new GT
		bool ** g_t = build_2D_array<bool>(vertexNum, false);
		int ** g_t_int = build_2D_array<int>(vertexNum, 0);

		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
			fscanf(fin, "%d", &g_t_int[i][j]);
			}
		}
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				g_t[i][j] = g_t_int[i][j] == 1 ? true : false;
			}
		}

#if 0
		printMatrix(vertexNum, g_t);
#endif

		fclose(fin);

		string tem_name = generateFileName(resutsDir, vertexNum, d, IterationNum,
			distribuion_type, SA_flag, workerNum, ratio ,  
			"-WQ-" + s_Worker_quality + "-Smooth-" + s_smothingType + "-rt-" + s_tasks_ratio + s_task_name, cur_time, ".txt");
		std::ofstream  of2(tem_name.c_str(), std::ofstream::out | std::ofstream::app);


		if (!of2.is_open())
			cout << "Not Open file " << tem_name << ".\n";

		vector <int> TruthIterationTimes(ExperimentTimes, 0);
		double
			SATD_tau = 0, // SA with truth discovery tau distance;
			CBT_tau = 0, // CBT tau distance;
			RC_tau = 0, // RC tau distance;
			GM_tau = 0, // GM tau distance;
			SATD_sec = 0, // SA with truth discovery running time;
			CBT_sec = 0, // CBT running time;
			RC_sec = 0, // RC running time;
			GM_sec = 0, // GM running time;
			SATD_WWM_tau = 0, // SATD with Worker Weight Modification or smoothing;
			SATD_WWM_sec = 0, // SATD with Worker Weight Modification or smoothing;
			// == true, using the TrustTransitivity, equation (7).
			// const bool IsTrustTransitivityFlag = true; 
			SATD_WWM_TD_sec = 0, // SATD_WWM running time, for truth discovery.
			SATD_WWM_Smoothing_sec = 0, // SATD_WWM running time, for smoothing.
			SATD_WWM_TC_sec = 0, // SATD_WWM running time, for Transitive closure.
			SATD_WWM_HP_sec = 0; // SATD_WWM running time, for finding best HP.

#if 1		  
			SA_SATD_CBT_RC_GM_Comparison(SATD_tau, CBT_tau, RC_tau, GM_tau, 
				SATD_sec, CBT_sec, RC_sec, GM_sec, SATD_WWM_tau, 
				SATD_WWM_sec, SATD_WWM_TD_sec, SATD_WWM_Smoothing_sec ,
				SATD_WWM_TC_sec, SATD_WWM_HP_sec,
				vertexNum, d, g_t, /*task assignment graph.*/
			ratio, coolrate, /*e.g. = (float) 0.95;*/
			IterationNum, /*e.g., = 100000;*/
			temperature, /*temperature, e.g., = 100000;*/
			ExperimentTimes, CBT_ExperimentTimes, permuStartingIdx, 
			distribuion_type, Worker_quality, smothingType,
			SA_flag, workerNum, tasks_ratio, prob_threshold, 
			Is_Beta_Simu_data_type, TruthIterationTimes, of2);
#endif

#if 1
			// Crowd_BT only
			/*w = 200; mean = 0; stddev = 0.001;*/
			int num_repeat = 1;
			int w = vertexNum *d *(tasks_ratio > 0 ? tasks_ratio : 1.0)
				* workerNum / (2 * num_repeat);
			if (w > 4000)
				w = 4000;
			int rep = 0;
			int * GTr_inver = new int[vertexNum];

			of2 << "\nGround truth For Crowd_BT is : \n";
			for (int i = 0; i < vertexNum; i++){
				GTr_inver[i] = vertexNum - 1 - i;
				of2 << i << " ";
			}

			of2 << "\n";
			vector<double> cbt_kt_gt(ExperimentTimes, 0.0);
			vector<double> cbt_elapse_time(ExperimentTimes, 0.0);
			if (rep < CBT_ExperimentTimes){
				// Active learning Crowd_BT
				cout << "Running Crowd_BT, workers = " << w << "," << rep + 1 <<
					"/" << CBT_ExperimentTimes << endl;
				of2 << "Running Crowd_BT, workers = " << w << "," << rep + 1 <<
					"/" << CBT_ExperimentTimes << endl;
				run_CrowdBT_for_SA_Comparison(
					cbt_kt_gt[rep], cbt_elapse_time[rep],
					vertexNum,
					// GTr, /*groud truth ranking;*/
					GTr_inver, /*reverse ground truth ranking, since we need ascending order for Crowd_BT.*/
					distribuion_type, Worker_quality, w, /*number of workers;*/
					num_repeat, Is_Beta_Simu_data_type);
				rep ++;
			}
			 
			CBT_sec = 0;
			// of1 << "\n CBT running running time in seconds: \n";
			for (int i = 0; i < CBT_ExperimentTimes; ++i){
				/*of1 << std::fixed << std::setprecision(8)
				<< v_cbt_elapse_time[i] << ", ";*/
				CBT_sec += cbt_elapse_time[i];
			}
			CBT_sec /= CBT_ExperimentTimes; // CBT running time;
			of2 << "\nCBT average running time is " << std::fixed << std::setprecision(8)
				<< CBT_sec << " seconds.\n";

			CBT_tau = 0;
			for (int i = 0; i < ExperimentTimes; ++i){
				CBT_tau += cbt_kt_gt[i];
			}

			CBT_tau /= CBT_ExperimentTimes; // CBT running time;
			of2 << " \nAverage KanTau(CBT, GroundTruth) = "
				<< std::fixed << std::setprecision(8) << CBT_tau << endl;
#endif

		// save the task assignment graph.

			of1 << "# " << gtfilelist[fidx] << std::fixed << std::setprecision(3) <<
				" truth discovery execution times are : ";
			for (int i = 0; i < ExperimentTimes; ++i){
				of1 << TruthIterationTimes[i] << ", ";
			} of1 << "\n";

			of1 << "# " << gtfilelist[fidx] << std::fixed << std::setprecision(3) <<
				" smoothing edges number is : ";
			for (int i = 0; i < ExperimentTimes; ++i){
				of1 << v_smoothing_counts[i] << ", allocation = " << v_smoothing_time_debug[0]
					<< ", process = " << v_smoothing_time_debug[1]
					<< ", small prob = " << v_smoothing_time_debug[2] << ";   ";
			} 

			of1 << "\n";

			of1 << std::fixed << std::setprecision(3)
				<< ratio << "\t" << std::fixed << std::setprecision(8)
				// accuracy
				<< SATD_tau << "\t"
				<< SATD_WWM_tau << "\t"
				<< RC_tau << "\t"
				<< GM_tau << "\t"
				<< CBT_tau << "\t"
				// time
				<< SATD_sec << "\t"
				<< SATD_WWM_sec << "\t"
				<< RC_sec << "\t"
				<< GM_sec << "\t"
				<< CBT_sec << "\t"
				<< vertexNum << "\t"
				<< SATD_WWM_TD_sec << "\t"
				<< SATD_WWM_Smoothing_sec << "\t"
				<< SATD_WWM_TC_sec << "\t"
				<< SATD_WWM_HP_sec << "\t"
				<< SATD_WWM_Smoothing_sec - v_smoothing_time_debug[0] << "\n";
				

	      // release memory;
		vector <int>().swap(TruthIterationTimes);
		delete_2D_Array<bool>(vertexNum, g_t);
		delete_2D_Array<int>(vertexNum, g_t_int);
	} /*end of each GT file*/
	of1.close();
	vector<double>().swap(v_tc_secs); 
	vector<double>().swap(v_td_secs); 
	vector<double>().swap(v_tc_secs); 
	return 0;
}
#endif
#endif