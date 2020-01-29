/**
* @file: Experiment2.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#define _CRT_SECURE_NO_WARNINGS
#include "Experiment.hpp"
#include "Baseline2.hpp"
#include "CrowdBT.hpp"
#include "truth_discovery.hpp"


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
	const int & workerNum,
	const int & SA_flag,
	const string & TaskName, // e.g., = SA, or  = SATD;
	std::ofstream & of1
	){

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" + to_string(now->tm_mday)
		+ "-" + to_string(now->tm_hour) + "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);
	
	// final experiment
	bool HPorNot = true;


	// display some important input parameters of main function
	std::cout << " \n//***************************" << endl;
	std::cout << " //*** Run " << TaskName << " From TC *********" << endl;
	std::cout << " //****************************" << endl;

	// tau distance;
	int * SAr = new int[vertexNum];

	TSPalgorithm alg(vertexNum, p_Gp_tc, permuStartingIdx);
	alg.Initialize(IterationNum, temperature, coolrate);
	// run SA;
	double score = 0;
	bool  ResultIsHP, iniIsHP;
	double iniScore = -1;
	alg.setMatrix(vertexNum, p_Gp_tc);
	// release memory;
	auto sa_begin = std::chrono::high_resolution_clock::now();
	std::vector<int> v_temp = alg.Run(SA_flag, 0, score, ResultIsHP,
		iniScore, iniIsHP, v_M);

	//printInitalTourFinalTour(iniScore, score, iniIsHP, ResultIsHP);
	alg.opt_path = make_pair(v_temp, score);
	for (auto i = 1; i < vertexNum; ++i){
		v_temp = alg.Run(SA_flag, i, score, ResultIsHP,
			iniScore, iniIsHP, v_M);

		//printInitalTourFinalTour(iniScore, score, iniIsHP, ResultIsHP);

		if (score < alg.opt_path.second)
			alg.opt_path = make_pair(v_temp, score);
	}
	auto sa_end = std::chrono::high_resolution_clock::now();
	// save the running time for SA;
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(sa_end - sa_begin).count();
	sec = dt;
	//alg.printHP(of1);
#if 1
	of1 << "\n" << TaskName << " result is : \n";
	std::cout << "\n" << TaskName << " result is : \n";
	for (auto i = 0; i < vertexNum - 1; ++ i){
		int v1 = alg.opt_path.first[i];
		int v2 = alg.opt_path.first[i+1];
		std::cout << v1 << " --( " << alg.tour->matrix[v1][v2] << " )-- ";
		of1 << v1 << " --( "       << alg.tour->matrix[v1][v2] << " )-- ";
	}
	std::cout << alg.opt_path.first[vertexNum - 1] << endl;
	of1 << alg.opt_path.first[vertexNum - 1] << endl;

#endif

	for (int i = 0; i < vertexNum; ++i)
		SAr[i] = alg.opt_path.first[i];

	std::cout << TaskName << " path score = " << alg.opt_path.second << endl;
	of1 << TaskName << " path score = " << alg.opt_path.second << endl;
	HPorNot = isHP<double, int>(p_Gp_tc, SAr, vertexNum);
	tau = kendallSmallN(SAr, TaskName/*e.g., = "SA"*/, GTr,
		"GroudTruth", vertexNum, false);


	of1 << "\n IsHP for each trial: ";
	if (HPorNot)
		of1 << "True, ";
	else
		of1 << "False, ";

	of1 << std::endl;

	of1 << "\n KanTau(" << TaskName << ", GroundTruth) = " << tau << endl;

	//***********************//
	// SA running time
	//***********************//
	//  running time;
	of1 << "\n" << TaskName << " running time = " << sec << " seconds.\n";
	//release 
	delete[] SAr;
}



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
	){
	auto t2 = std::chrono::high_resolution_clock::now();
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" + to_string(now->tm_mday)
		+ "-" + to_string(now->tm_hour) + "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);

	// display some important input parameters of main function
	std::cout << " \n//************************************" << endl;
	std::cout << " //*** Run GM From Voting Matrix *******" << endl;
	std::cout << " //*************************************" << endl;

	string title = "Running GM From Voting Matrix";
	saveParameters(title, vertexNum, d, distribuion_type, w,
		cur_time, of1);

	// tau distance;
	int * GMr = new int[vertexNum];

	
	
	for (int i = 0; i < 5000; ++i){
		//of1 << "\nGround truth is : \n";
		for (int i = 0; i < vertexNum; i++){
			GMr[i] = i;
		}
		CondorcetFuseSort(GMr, 0, vertexNum - 1, vertexNum,
		w, // each of the w search systems, or w voters; 
		matrice // each voter's voting result;
		);
	}
	

	tau = kendallSmallN(GMr, "GM", GTr,
		"GroundTruth", vertexNum, false);
	auto t3 = std::chrono::high_resolution_clock::now();
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
	sec = dt;
	std::cout << "Finish GM in " << sec << " seconds" << endl;
	//of1 << "Finish RC in " << dt << " seconds" << endl;


	of1 << "\n KanTau(GM, GT)  = " << tau << "\n";
	of1 << "\n GM running running time = " << sec << " seconds.\n";
	//release 
	delete[] GMr;
}

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
	){
	auto t2 = std::chrono::high_resolution_clock::now();
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" +
		to_string(now->tm_mday) + "-" + to_string(now->tm_hour)
		+ "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);


	// display some important input parameters of main function
	std::cout << " \n//************************************" << endl;
	std::cout << " //*** Run RC  From Voting Matrix ******" << endl;
	std::cout << " //*************************************" << endl;
	string title = "Running RC";
	saveParameters(title, vertexNum, d, distribuion_type, w,
		cur_time, of1);

	// tau distance;
	int * RCr = new int[vertexNum];

		//run RC and accumulate dis
		RepeatChoice rc(vertexNum, matrice);
		
		vector<int> ranking_RC = rc.process();
		
		for (int i = 0; i < vertexNum; ++i){
			RCr[i] = ranking_RC[i];
		}

		tau = kendallSmallN(RCr, "RC", GTr,
			"GroundTruth", vertexNum, false);
		auto t3 = std::chrono::high_resolution_clock::now();
		auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
		sec = dt;
		cout << "Finish RC in " << sec << " seconds" << endl;
		//of1 << "Finish RC in " << dt << " seconds" << endl;
		vector<int>().swap(ranking_RC);
	
	of1 << "KanTau(RC, GroundTruth) =  " << tau << endl;
	of1 << "RC running running time = " << sec <<" seconds.\n";
	//release memory
	delete[] RCr;
}


