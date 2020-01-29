/**
* @file: Experiment.cpp
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

// For task 2, compare TA and SA, in terms of running time and 
// accuracy (that is Kendall tau distance.)

// For task 2, compare TA and SA, in terms of running time and 
// accuracy (that is Kendall tau distance.)

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
	){

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" + to_string(now->tm_mday)
		+ "-" + to_string(now->tm_hour) + "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);

	int w = 20; // worker number;
	int e = vertexNum * d / 2; // edge number;

	// Kendall Tau distance (SA, Ground_Truth);
	vector<double> ta_kt_gt_reverse(ExperimentTimes, 0.0);
	vector<double> v_ta_elapse_time(ExperimentTimes, 0.0);

	// display some important input parameters of main function
	std::cout << " \n//*****************************************" << endl;
	std::cout << " //*** TA and Ground Truth Comparison! *****" << endl;
	std::cout << " //*****************************************" << endl;

	// display some important input parameters of main function
	of1 << " //*****************************************" << endl;
	of1 << " //*** TA and Ground Truth Comparison! *****" << endl;
	of1 << " //*****************************************" << endl;

	// save those parameters
	of1
		<< "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl
		<< "Experiment Times for average result = " << ExperimentTimes << endl;
	if (distribuion_type == 0)
		of1 << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		of1 << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	of1 << "Parallel vertex index if needed = " << parallelIdx << endl
		<< "Current time is " << cur_time << endl;

	// tau distance;
	int * TAr = new int[vertexNum];
	int * GTr_inver = new int[vertexNum]; // here we use inverse, since 
	// they just use opposite "preference" direction.

	double ** p_Gp_tc = new double*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		p_Gp_tc[i] = new double[vertexNum];
	}

	std::map<int, int> invert_index;
	vector<int> ground_truth, vertexSet;
	for (int i = 0; i < vertexNum; i++){
		ground_truth.push_back(i);
		vertexSet.push_back(i);
	}
	vector<double > v_weights(vertexNum*vertexNum, 0.0);
	/* initialize random seed: */
	/* initialize random seed: */
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);

	for (int rep = 0; rep < ExperimentTimes; ++rep){ // multiple experiment trials;
		clock_t be1 = clock();
		// randomize the ground truth;
		std::shuffle(ground_truth.begin(), ground_truth.end(), engine);

		of1 << "\nGround truth is : \n";
		for (int i = 0; i < vertexNum; i++){
			invert_index[ground_truth[i]] = i;
			GTr_inver[i] = ground_truth[vertexNum - 1 - i];
			of1 << GTr_inver[i] << " ";
		}
		of1 << "\n";
		// generate the normalized preference graph G_P;
		// that is, G_P(i, j) + G_P(j,i) = 1, where i != j;
		// G_P(i, i) = 0;
		std::vector<std::vector<double>> gp =
			generate_aggregated_matrix(vertexNum, w, g_t,
			invert_index, distribuion_type, 
			quality, // different levels of error-rate for the worker's quality.
			max_per_worker,
			smothingType);

		// using G_P to generate the transitive closure G_P*;
		std::vector<std::vector<double>> gp_tc = generateTransitionMatrix_Eigen(gp);


		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				v_weights[i*vertexNum + j] = gp_tc[i][j];
		// run TA;
		// initialize static variables in class TA;
		clock_t ta_begin = clock();
		TA ta(vertexNum, e, 5 /*top_k*/, false, false);
		ta.initiaWeights(v_weights);
		//ta.getSortedWeights(v_weights);
		// if NO-HP;
		of1 << "TA result is :\n";
		if (!ta.bruteForce(vertexSet, of1)){
			rep--;
			// if no-HP,
			// directly continue next experiment trial;
			// without doing the following SA at all;
			continue;
		}
		else{

			TAr = ta.getmaxPath();
			clock_t ta_end = clock();
			// save the running time for TA
			v_ta_elapse_time[rep] = double(ta_end - ta_begin) / CLOCKS_PER_SEC;
			for (int i = 0; i < vertexNum; i++){
				std::cout << TAr[i] << ", ";
				of1 << TAr[i] << ", ";
			}
			ta_kt_gt_reverse[rep] = kendallSmallN(TAr, "TA", GTr_inver,
				"Gt_inver", vertexNum, false);
		}

		// release memory;
		std::vector<std::vector<double>>().swap(gp);
		std::vector<std::vector<double>>().swap(gp_tc);
	}/*end of each time experiment*/

	//***********************//
	// accuracy
	//***********************//
	// TA: calculate tau-distance on average;
	std::sort(ta_kt_gt_reverse.begin(), ta_kt_gt_reverse.end());
	double ta_sco = .0;

	if (ExperimentTimes > 2){
		for (int i = 1; i < ExperimentTimes - 1; ++i){
			ta_sco += ta_kt_gt_reverse[i];
		}
		ta_sco /= (ExperimentTimes - 2);
	}

	else{
		for (int i = 0; i < ExperimentTimes; ++i){
			ta_sco += ta_kt_gt_reverse[i];
		}

		ta_sco /= (ExperimentTimes);
	}

	of1 << "\n KanTau(TA, GT_RES) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << ta_kt_gt_reverse[i] << ", ";
	}
	of1 << " \nAverage KanTau(TA, GT_RES) = " << ta_sco << endl;


	//***********************//
	// running time
	//***********************//
	// TA running time;
	double ta_time = 0, sa_time = 0;
	of1 << "\n TA running time in seconds: \n";
	for (int i = 0; i < v_ta_elapse_time.size(); ++i){
		of1 << v_ta_elapse_time[i] << ", ";
		ta_time += v_ta_elapse_time[i];
	}
	ta_time /= ExperimentTimes;
	of1 << "\nTA average time is " << ta_time << " seconds.\n";

	//of1.close();
	//release 
	delete[] TAr;
	delete[] GTr_inver;

	for (int i = 0; i < vertexNum; ++i){
		delete[] p_Gp_tc[i];
	}
	delete p_Gp_tc;

	vector<double >().swap(v_weights);
	vector<int>().swap(ground_truth);
	vector<int>().swap(vertexSet);
	vector<double>().swap(ta_kt_gt_reverse);
	vector<double>().swap(v_ta_elapse_time);
}

void set_GT_False(const int& n, bool **GT){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			GT[i][j] = false;
}

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
	){
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" + to_string(now->tm_mday)
		+ "-" + to_string(now->tm_hour) + "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);

	// final experiment
	// G_T generation;
	int w = 20; // worker number;

	// Kendall Tau distance (SA, Ground_Truth);
	vector<double> sa_kt_ta(ExperimentTimes, 0.0); // SA and TA
	vector<double> sa_kt_gt_reverse(ExperimentTimes, 0.0);
	vector<double> ta_kt_gt_reverse(ExperimentTimes, 0.0);
	vector<double> v_ta_elapse_time(ExperimentTimes, 0.0);
	vector<double> v_sa_elapse_time(ExperimentTimes, 0.0);


	// display some important input parameters of main function
	std::cout << " \n//***********************************" << endl;
	std::cout << " //*** TA and SA Comparison! *********" << endl;
	std::cout << " //***********************************" << endl;

	// display some important input parameters of main function
	of1 << " \n//***********************************" << endl;
	of1 << " //*** TA and SA Comparison! *********" << endl;
	of1 << " //***********************************" << endl;

	// save those parameters
	of1 << "vertexNum = " << vertexNum << endl
		<< "vertex degree = " << d << endl
		<< "Initial Temperature T = " << temperature << endl
		<< "Cooling rate = " << coolrate << endl
		<< "Iteration Number = " << IterationNum << endl
		<< "Experiment Times for average result = " << ExperimentTimes << endl;
	if (distribuion_type == 0)
		of1 << "Simulated Data Distribution type = UNIFORM_DISTRIBUTION\n";
	else if (distribuion_type == 1)
		of1 << "Simulated Data Distribution type = NORMAL_DISTRIBUTION\n";
	of1 << "Parallel vertex index if needed = " << parallelIdx << endl
		<< "Current time is " << cur_time << endl;


	// tau distance;
	int * TAr = new int[vertexNum];
	int * SAr = new int[vertexNum];
	int * GTr_inver = new int[vertexNum]; // here we use inverse, since 
	// they just use opposite "preference" direction.

	double ** p_Gp_tc = new double*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		p_Gp_tc[i] = new double[vertexNum];
	}

	std::map<int, int> invert_index;
	vector<int> ground_truth, vertexSet;
	for (int i = 0; i < vertexNum; i++){
		ground_truth.push_back(i);
		vertexSet.push_back(i);
	}

	int e = vertexNum * d / 2;
	/* initialize random seed: */
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	auto engine = std::default_random_engine(seed);

	vector<double > v_weights(vertexNum*vertexNum, 0.0);
	vector<vector<vector<int> >> matrice;
	std::vector<std::vector<double>> gp;
	// using G_P to generate the transitive closure G_P*;
	std::vector<std::vector<double>> gp_tc;

	// multiple experiment trials;
	for (int rep = 0; rep < ExperimentTimes; ++rep){
		clock_t be1 = clock();
		// randomize the ground truth;
		std::shuffle(ground_truth.begin(), ground_truth.end(), engine);

		of1 << "\nGround truth is : \n";
		for (int i = 0; i < vertexNum; i++){
			invert_index[ground_truth[i]] = i;
			GTr_inver[i] = ground_truth[vertexNum - 1 - i];
			of1 << GTr_inver[i] << " ";
		}
		of1 << "\n";
		
		matrice = generate_sim_votes(vertexNum, workerNum, g_t, invert_index,
			distribuion_type,
			quality, // different levels of error-rate for the worker's quality.
			max_per_worker);

		gp =generate_aggregated_matrix(matrice, smothingType);
		// using G_P to generate the transitive closure G_P*;
		gp_tc = generateTransitionMatrix_Eigen(gp);

		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				v_weights[i*vertexNum + j] = gp_tc[i][j];
		// run TA;
		// initialize static variables in class TA;
		clock_t ta_begin = clock();
		TA ta(vertexNum, e, 5 /*top_k*/, false, false);
		ta.initiaWeights(v_weights);
		//ta.getSortedWeights(v_weights);
		// if NO-HP;
		of1 << "TA result is :\n";
		if (!ta.bruteForce(vertexSet, of1)){
			rep--;
			// if no-HP,
			// directly continue next experiment trial;
			// without doing the following SA at all;
			continue;
		}
		else{

			TAr = ta.getmaxPath();
			clock_t ta_end = clock();
			// save the running time for TA
			v_ta_elapse_time[rep] = double(ta_end - ta_begin)
				/ CLOCKS_PER_SEC;
			for (int i = 0; i < vertexNum; i++){
				std::cout << TAr[i] << ", ";
				of1 << TAr[i] << ", ";
			}
			ta_kt_gt_reverse[rep] = kendallSmallN(TAr, "TA",
				GTr_inver, "Gt_inver", vertexNum, false);
		}

		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				p_Gp_tc[i][j] = gp_tc[i][j];
			}
		}

		TSPalgorithm alg(vertexNum, p_Gp_tc, permuStartingIdx);
		alg.Initialize(IterationNum, temperature, coolrate);
		double score = 0;
		clock_t sa_begin = clock();
		
		std::vector<int> v_temp = alg.Run(SA_flag, score, 0, matrice);
		alg.opt_path = make_pair(v_temp, score);
		for (auto i = 1; i < vertexNum; ++i){
			v_temp = alg.Run(SA_flag, score, i, matrice);
			if (score < alg.opt_path.second)
				alg.opt_path = make_pair(v_temp, score);
		}
		clock_t sa_end = clock();

		// save the running time for SA;
		v_sa_elapse_time[rep] = double(sa_end - sa_begin) / CLOCKS_PER_SEC;
		// alg.printHP(of1);
		of1 << "\nSA result is :\n";
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
		std::cout << "\nSA result is :\n";
#endif
		for (auto i : alg.opt_path.first){
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
			std::cout << i << ", ";
#endif
			of1 << i << ", ";
		}
		for (int i = 0; i < vertexNum; ++i)
			SAr[i] = alg.opt_path.first[i];

		std::cout << "SA path score = " << alg.opt_path.second << endl;
		of1 << "SA path score = " << alg.opt_path.second << endl;
		sa_kt_gt_reverse[rep] = kendallSmallN(SAr, "SA", GTr_inver,
			"Gt_inver", vertexNum, false);
		sa_kt_ta[rep] = kendallSmallN(SAr, "SA", TAr, "TA", vertexNum, false);

	}/*end of each time experiment*/


	// release memory;
	vector<vector<vector<int> >>().swap(matrice);
	std::vector<std::vector<double>>().swap(gp);
	std::vector<std::vector<double>>().swap(gp_tc);

	//***********************//
	// accuracy
	//***********************//
	// TA: calculate tau-distance on average;
	std::sort(ta_kt_gt_reverse.begin(), ta_kt_gt_reverse.end());
	// SA: calculate tau-distance on average;
	std::sort(sa_kt_ta.begin(), sa_kt_ta.end());
	std::sort(sa_kt_gt_reverse.begin(), sa_kt_gt_reverse.end());
	double ta_sco = .0, sa_sco = 0.0, ta_sa_sco = 0.0;

	if (ExperimentTimes > 2){
		for (int i = 1; i < ExperimentTimes - 1; ++i){
			sa_sco += sa_kt_gt_reverse[i];
			ta_sco += ta_kt_gt_reverse[i];
			ta_sa_sco += sa_kt_ta[i];
		}
		sa_sco /= (ExperimentTimes - 2);
		ta_sco /= (ExperimentTimes - 2);
		ta_sa_sco /= (ExperimentTimes - 2);
	}

	else{
		for (int i = 0; i < ExperimentTimes; ++i){
			sa_sco += sa_kt_gt_reverse[i];
			ta_sco += ta_kt_gt_reverse[i];
			ta_sa_sco += sa_kt_ta[i];
		}
		sa_sco /= (ExperimentTimes);
		ta_sco /= (ExperimentTimes);
		ta_sa_sco /= (ExperimentTimes);
	}

	of1 << "\n KanTau(TA, GT_RES) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << ta_kt_gt_reverse[i] << ", ";
	}
	of1 << " \nAverage KanTau(TA, GT_RES) = " << ta_sco << endl;


	of1 << "\n KanTau(SA, GT_RES) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << sa_kt_gt_reverse[i] << ", ";
	}
	of1 << " \nAverage KanTau(SA, GT_RES) = " << sa_sco << endl;

	of1 << "\n KanTau(SA, TA) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << sa_kt_ta[i] << ", ";
	}
	of1 << " \nAverage KanTau(SA, TA) = " << ta_sa_sco << endl;

	//***********************//
	// running time
	//***********************//
	// TA running time;
	double ta_time = 0, sa_time = 0;
	of1 << "\n TA running time in seconds: \n";
	for (int i = 0; i < v_ta_elapse_time.size(); ++i){
		of1 << v_ta_elapse_time[i] << ", ";
		ta_time += v_ta_elapse_time[i];
	}
	ta_time /= ExperimentTimes;
	of1 << "\nTA average time is " << ta_time << " seconds.\n";

	of1 << "\n SA running running time in seconds: \n";
	for (int i = 0; i < v_sa_elapse_time.size(); ++i){
		of1 << v_sa_elapse_time[i] << ", ";
		sa_time += v_sa_elapse_time[i];
	}
	sa_time /= ExperimentTimes;
	of1 << "\nSA average running time is " << sa_time << " seconds.\n";

	//of1.close();
	//release 
	delete[] TAr;
	delete[] SAr;
	delete[] GTr_inver;

	for (int i = 0; i < vertexNum; ++i){
		delete[] p_Gp_tc[i];
	}
	delete p_Gp_tc;

	vector<double >().swap(v_weights);
	vector<int>().swap(ground_truth);
	vector<int>().swap(vertexSet);
	vector<double>().swap(sa_kt_ta);
	vector<double>().swap(sa_kt_gt_reverse);
	vector<double>().swap(ta_kt_gt_reverse);
	vector<double>().swap(v_ta_elapse_time);
	vector<double>().swap(v_sa_elapse_time);
}

// try to do task assignment to generate G_T graph;
// if the degree d is not appropriate,
// we will try some other degree within some time threshold;
// otherwise, false will be returned.
bool task_assignment_graph(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	int & d// degree 
	){

	// task assignment G_T;
	vector<int> v_degrees;
	v_degrees.push_back(d);
	int max_times = TASK_ASSIGNMENT_MAX_TIMES;
	for (int i = 1; i < max_times; ++i){
		if (d + 2 * i < vertexNum - 1)
			v_degrees.push_back(d + 2 * i);
		if (d - 2 * i > 2)
			v_degrees.push_back(d - 2 * i);
	}

	double thre = THRESHOLD_TIME_IN_SECONDS;
	bool flg = false;
	for (auto temp_d : v_degrees){
		// task assignment succeed;
		if (flg = get_G_T(g_t, vertexNum, thre, temp_d)){
			d = temp_d;
			cout << " Success! Task assignment G_T for degree d = "
				<< temp_d << "\n";
			/*of1 << " Success! Task assignment G_T for degree d = "
				<< temp_d << "\n";*/
			break;
		}
	}

	if (!flg){
		cout << "Error! Task assignment failed! No TA done!\n";
		/*of1 << "Error! Task assignment failed! No TA done!\n";*/
		return false;
	}
	else
		return true;
}

bool get_G_T(
	bool ** g_t, // undirected task assignment graph G_T;
	const int & vertexNum, // number of vertex ;
	const double & thre, // time threshold, since we do not
	// want the task assignment process falls into 
	// an endless loop.
	const int & d// degree 
	){
	// since sometimes some degree d cannot 
	// guarantee the task assignment is possible.
	// thus, we set a time threshold to test whether 
	// the task assignment can be done.
	double delta_time = .0; // duration in seconds;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	// 1.e-9: nano-seconds to seconds;
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

	bool isDirect = false; // is directed graph or not;
	int pathNodeLength = 4; // it does not matter for this experiment;
	int CircleNodeLength = 4;// it does not matter for this experiment; 
	myGraph gt(vertexNum, isDirect, d, vertexNum - 1,
		pathNodeLength, CircleNodeLength);
	// reset the matrix elements to falses(i.e., zeros);
	gt.resetMatrix();
	bool re = true;
	while (!(re = gt.taskAssign(thre))){
		// reset the matrix elements to falses(i.e., zeros);
		cout << "Reset the matrix elements to zeros,"
			<< " and degrees to zeros,  for next round G_T generation.\n";
		gt.initialize();
		t2 = std::chrono::high_resolution_clock::now();
		// 1.e-9: nano-seconds to seconds;
		dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
		if (dt > thre){
			std::cout << "Time too long, break the loop!\n";
			break;
		}
	}

	if (re)
		//copy G_T
		gt.copyG_T(g_t);
	else
		g_t = NULL;
	return re;
}



std::vector<std::vector<double>> getGTCfromGPInFile(const string path,
	const int & vertexNum){
	std::vector<CSVRow> v_gp = readCSVFiles(path);
	if (v_gp.empty()){
		cout << "Failed in reading GP matrix from" << path << ".\n";
		return std::vector<std::vector<double>>(0, vector<double>(0,0));
	}
	else
	{
		std::vector<std::vector<double>> gp(vertexNum, vector<double>(vertexNum, .0));
		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				gp[i][j] = stod(v_gp[0][i * vertexNum + j]);
		std::vector<std::vector<double>> gtc = generateTransitionMatrix_Eigen(gp);
		return gtc;
	}
}

bool getGTCfromGPInFile(const string path, double ** p_gtc,
	const int & vertexNum){

	std::vector<CSVRow> v_gp = readCSVFiles(path);
	if (v_gp.empty()){
		cout << "Failed in reading GP matrix from" << path << ".\n";
		return false;
	}
	else
	{
		std::vector<std::vector<double>> gp(vertexNum, vector<double>(vertexNum, .0));
		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				gp[i][j] = stod(v_gp[i][j]);
		std::vector<std::vector<double>> gtc = generateTransitionMatrix_Eigen(gp);
		for (int i = 0; i < vertexNum; ++i)
			for (int j = 0; j < vertexNum; ++j)
				p_gtc[i][j] = gtc[i][j];
		// release memory
		std::vector<CSVRow>().swap(v_gp);
		std::vector<std::vector<double>>().swap(gp);
		std::vector<std::vector<double>>().swap(gtc);
		return true;
	}
}


bool Run_Make_GT(
	const int & vertexNum,
	int & d,
	bool ** g_t,
	ofstream & of
	){
	std::cout << " //**********************" << endl;
	std::cout << " //*** Generate G_T *****" << endl;
	std::cout << " //**********************" << endl;
	// set zero values;
	set_GT_False(vertexNum, g_t);
	// if not successful;
	if (!task_assignment_graph(g_t, vertexNum, d)){
		std::cout << " Error! Task Assignment failed "
			<< " for v = " << vertexNum
			<< " , d = " << d << endl;
		of << " Error! Task Assignment failed "
			<< " for v = " << vertexNum
			<< " , d = " << d << endl;
		return false;
	}
	else
		return true;
}