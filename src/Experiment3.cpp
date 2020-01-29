/**
* @file: Experiment3.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#define _CRT_SECURE_NO_WARNINGS
#include "Experiment.hpp"
#include "Baseline2.hpp"
#include <iomanip> /*IO Manipulators*/
#include "CrowdBT.hpp"
#include "truth_discovery.hpp"

/*a threshold of the small_prob to be generated.*/
extern double prob_threshold;
extern double _ACCURCY;
extern int  totalTime;
extern vector<double> v_tc_secs;
extern vector<double> v_td_secs; // truth discovery;
extern vector<int> v_smoothing_counts;
extern int temp_main;
extern vector<double> v_smoothing_time_debug;

// compare all the algorithms using the same simulation data.
// SA2: our algorithm, with weights updaing used in truth discovery model;
// CBT: Crowd_BT modle;
// RC: Repeated Choice;
// GM: Condorcet Fuse Sort, or Quick Sort (QS);
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
	double & SATD_WWM_sec, // SATD_WWM running time;
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
	){

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-" + to_string(now->tm_mday)
		+ "-" + to_string(now->tm_hour) + "-" + to_string(now->tm_min)
		+ "-" + to_string(now->tm_sec);

	// Kendall Tau distance (SA, Ground_Truth);
	vector<double> satd_kt_gt(ExperimentTimes, 0.0);
	vector<double> cbt_kt_gt(ExperimentTimes, 0.0);
	vector<double> rc_kt_gt(ExperimentTimes, 0.0);
	vector<double> gm_kt_gt(ExperimentTimes, 0.0);
	vector<double> satd_elapse_time(ExperimentTimes, 0.0);
	vector<double> cbt_elapse_time(ExperimentTimes, 0.0);
	vector<double> rc_elapse_time(ExperimentTimes, 0.0);
	vector<double> gm_elapse_time(ExperimentTimes, 0.0);
	vector<double> satd_wwm_elapse_time(ExperimentTimes, 0.0);
	vector<double> satd_wwm_kt_gt(ExperimentTimes, 0.0);
	//vector <int> TruthIterationTimes(ExperimentTimes, 0);
	

	// display some important input parameters of main function
	std::cout << "\n//*****************************************************" << endl;
	std::cout <<   "//*** Compare SA, SATD, CBT, RC and GM Together! *****" << endl;
	std::cout <<   "//*****************************************************" << endl;

	// tau distance;
	int * GTr_inver = new int[vertexNum];
	int * GTr = new int[vertexNum];

	double ** p_Gp_tc = build_2D_array<double>(vertexNum, .0); 

	std::map<int, int> invert_index;
	vector<int> ground_truth;
	for (int i = 0; i < vertexNum; i++){
		ground_truth.push_back(i);
	}

	
	vector<vector<vector<int>>> matrice;
	//std::vector<std::vector<double>> gp;
	// using G_P to generate the transitive closure G_P*;
	// std::vector<std::vector<double>> gp_tc;

	std::vector<std::vector<double>> gp_SATD; 
	std::vector<std::vector<double>> gp_tc_SATD;
	// wwm: worker weight modification.
	std::vector<std::vector<double>> gp_SATD_wwm;
	std::vector<std::vector<double>> gp_tc_SATD_wwm;

	//*****************************
	// Truth discovery struct
	//*****************************
	/* In the current version, I defined all the arrays in the dynamic way 
	 * as you required. But another problem shows up, allocating and freeing 
	 * such huge amount of memory frequently will also spend lots of time.
	 * For this reason, the init function now has a parameter called "first_time", 
	 * it is a boolean variable (0 or 1). When it's 0, the program will 
	 * not allocate new memory but just clean it, and you can use the 
	 * existing arrays to do you consist experiment. If it's 1, 
	 * the program will allocate memory firstly and then clean it.
	 * So, when you start my program at the first time, you must ** execute init(1)**.
	 * And for the later exp(**10 time for the same size of data**), you **exec init(0)**.
	 */
	
	td_wieghts::truth_discovery * td = new td_wieghts::truth_discovery();
	int max_task_num = vertexNum * d / 2 + 1000;
	td->init(1, vertexNum, workerNum, max_task_num);

	/* initialize random seed: */
	double t_voting_matrix = 0.0;
	double t_satd_tc = 0.0;

	// Separate step running time for SATD_WWM.
	    SATD_WWM_TC_sec = 0, // SATD_WWM running time, for Transitive closure.
		SATD_WWM_TD_sec = 0, // SATD_WWM running time, for truth discovery.
		SATD_WWM_Smoothing_sec = 0, // SATD_WWM running time, for smoothing.
		SATD_WWM_HP_sec = 0; // SATD_WWM running time, for finding best HP.

	// multiple experiment trials;
	for (int rep = 0; rep < ExperimentTimes; ++rep){
		if (rep % 2 == 0)
			cout << "Repeat " << rep << "/" << ExperimentTimes
			<< " trails.\n";

		// randomize the ground truth;
		// std::random_shuffle(ground_truth.begin(), ground_truth.end());

		/*cout << "Ground Truth is : \n";
		for (auto i : ground_truth)
		cout << i << " > ";
		cout << endl;*/

		of1 << "\nGround truth is : \n";
		for (int i = 0; i < vertexNum; i++){
			invert_index[ground_truth[i]] = i;
			GTr[i] = ground_truth[i];
			GTr_inver[i] = ground_truth[vertexNum - 1 - i];
			of1 << GTr[i] << " ";
		}
		of1 << "\n";

		int edge_num = get_edge_num(g_t, vertexNum);
		// maximum number of pairwise comparison per worker can do.
		int max_per_worker = floor((double)edge_num * tasks_ratio);
		if (max_per_worker < 3) max_per_worker = 3;

		auto t_start = std::chrono::high_resolution_clock::now();

		matrice = generate_sim_votes(vertexNum, workerNum,
			g_t, invert_index, distribuion_type,quality, max_per_worker);

		//gp = generate_aggregated_matrix(matrice, smothingType);
		auto t_end = std::chrono::high_resolution_clock::now();
		t_voting_matrix += 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();


		// time counting;
		t_start = std::chrono::high_resolution_clock::now();
		if (!gp_SATD.empty())
			gp_SATD.clear();
		vector<double> worker_weights(workerNum, 0.0);

		
		vector<vector <td_wieghts::client_vote>> v_client_votes =
			generate_aggregated_matrix_SATD(g_t, matrice, gp_SATD,
			worker_weights, td, max_task_num, TruthIterationTimes[rep], smothingType);


#if 0
		// store the tasks;
		vector <pair<int, int>> v_tasks;
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < i; ++j){
				if (g_t[i][j] == true){
					v_tasks.push_back(make_pair(i, j));
				}
			}
		}
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

#endif
		
		t_end = std::chrono::high_resolution_clock::now();
		SATD_WWM_TD_sec += 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
		v_td_secs[rep] = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
		const double mean = 0.0;
		const double stddev_scalar = 1.0;
		// apply the smoothing to the preference graph.


		t_start = std::chrono::high_resolution_clock::now();
		gp_SATD_wwm = gp_SATD;
		Smoothing_UpdateProb_with_Workers_Quality(
			worker_weights, g_t, vertexNum, prob_threshold,
			mean, stddev_scalar, v_client_votes,
			gp_SATD_wwm);
		t_end = std::chrono::high_resolution_clock::now();
		SATD_WWM_Smoothing_sec += 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
		v_smoothing_counts[rep] = v_smoothing_time_debug[3];

		t_start = std::chrono::high_resolution_clock::now();
		gp_tc_SATD_wwm = generateTransitionMatrix_Eigen(gp_SATD_wwm);
		t_end = std::chrono::high_resolution_clock::now();
		v_tc_secs[rep] = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
		SATD_WWM_TC_sec += 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();


		// run RC (Repeated Choice)
		Run_RC_From_VotingMatrix(
			rc_kt_gt[rep], rc_elapse_time[rep],
			// GTr, // groud truth ranking;
			GTr_inver, // reverse groud truth ranking;
			matrice, vertexNum, d, ExperimentTimes, distribuion_type,
			workerNum, of1);

		//Run GM
		Run_GM_From_VotingMatrix(gm_kt_gt[rep], gm_elapse_time[rep],
			// GTr, // groud truth ranking;
			GTr_inver, // reverse groud truth ranking;
			matrice, vertexNum, d, distribuion_type, workerNum, of1);


#if 1
#if 0
		if (rep < CBT_ExperimentTimes){
			of1 << "GT is:\n";
			fprint_2D<bool>(vertexNum, g_t, of1);
			
			of1 << "v_client_votes (saved for Crowd_BT) is:\n";
			for (int i = 0, s = matrice.size(); i < s; ++i){
				of1 << " Worker " << i << " votes:\n";
				int s2 = v_client_votes[i].size();
				for (int j = 0; j < s2 - 1; ++j){
					of1 << v_client_votes[i][j].task_id << " "
						<< v_client_votes[i][j].who << " ";
				}
				of1 << v_client_votes[i][s2-1].task_id << " "
					<< v_client_votes[i][s2 -1].who << "\n";
			}
		}

		//of1 << "GP_SATD is:\n";
		//fprint_2D_vector<double>(vertexNum, gp_SATD, of1);
		//of1 << "GP_SATD_wwm (wwm: worker weight modification) is:\n";
		//fprint_2D_vector<double>(vertexNum, gp_SATD_wwm, of1);
		//of1 << "GTC_SATD_wwm (wwm: worker weight modification) is:\n";
		//fprint_2D_vector<double>(vertexNum, gp_tc_SATD_wwm, of1);
#endif
#if 0	
		//Run SATD From TC
		// run SATD, with truth discovery;
		t_start = std::chrono::high_resolution_clock::now();
		gp_tc_SATD = generateTransitionMatrix_Eigen(gp_SATD);
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				p_Gp_tc[i][j] = gp_tc_SATD[i][j];
			}
		}
		t_end = std::chrono::high_resolution_clock::now();
		t_satd_tc += 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();

		Run_SA_SATD_From_TC(satd_kt_gt[rep], satd_elapse_time[rep],
			// GTr_inver, // reverse groud truth ranking;
			GTr, // groud truth ranking;
			matrice,
			p_Gp_tc, // transitive closure matrix G_TC;
			vertexNum, d, coolrate, IterationNum, temperature,
			permuStartingIdx, distribuion_type,
			workerNum, SA_flag, "SATD", of1);


		if (rep < CBT_ExperimentTimes){
			of1 << "GT is:\n";
			fprint_2D<bool>(vertexNum, g_t, of1);
			of1 << "GP (saved for Crowd_BT) is:\n";
			fprint_2D_vector<double>(vertexNum, gp, of1);
		}

		of1 << "GP_SATD is:\n";
		fprint_2D_vector<double>(vertexNum, gp_SATD, of1);
		of1 << "GTC_SATD is:\n";
		fprint_2D_vector<double>(vertexNum, gp_tc_SATD, of1);

		of1 << "GP_SATD_wwm (wwm: worker weight modification) is:\n";
		fprint_2D_vector<double>(vertexNum, gp_SATD_wwm, of1);
		of1 << "GTC_SATD_wwm (wwm: worker weight modification) is:\n";
		fprint_2D_vector<double>(vertexNum, gp_tc_SATD_wwm, of1);
#endif



		//Run SATD_TT From TC
		for (int i = 0; i < vertexNum; ++i){
			for (int j = 0; j < vertexNum; ++j){
				p_Gp_tc[i][j] = gp_tc_SATD_wwm[i][j];
			}
		}
		Run_SA_SATD_From_TC(satd_wwm_kt_gt[rep], satd_wwm_elapse_time[rep],
			// GTr_inver, // reverse groud truth ranking;
			GTr, // groud truth ranking;
			matrice,
			p_Gp_tc, // transitive closure matrix G_TC;
			vertexNum, d, coolrate, IterationNum, temperature,
			permuStartingIdx, distribuion_type,
			workerNum, SA_flag, "SATD_WWM", of1);


#endif

#if 0
		/*w = 200; mean = 0; stddev = 0.001;*/
		int num_repeat = 4;
		int w = vertexNum *d *(tasks_ratio > 0 ? tasks_ratio : 1.0)
			* workerNum / (2 * num_repeat);
		if (w > 4000)
			w = 4000;
		if (rep < CBT_ExperimentTimes){
			// Active learning Crowd_BT
			cout << "Running Crowd_BT, workers = " << w << "," << rep + 1 <<
				"/" << CBT_ExperimentTimes << endl;
			of1 << "Running Crowd_BT, workers = " << w << "," << rep + 1 <<
				"/" << CBT_ExperimentTimes << endl;
			run_CrowdBT_for_SA_Comparison(
			cbt_kt_gt[rep], cbt_elapse_time[rep],
			vertexNum,
			// GTr, /*groud truth ranking;*/
			GTr_inver, /*reverse ground truth ranking, since we need ascending order for Crowd_BT.*/
			distribuion_type, quality, w, /*number of workers;*/
			num_repeat, Is_Beta_Simu_data_type);
		}
#endif

	}/*end of each time experiment*/

	// release memeory timely;
	delete_2D_Array<double>(vertexNum, p_Gp_tc);
	delete[] GTr;
	delete[] GTr_inver;
	vector<int>().swap(ground_truth);
	vector<vector<vector<int>>>().swap(matrice);
	//std::vector<std::vector<double>>().swap(gp);
	std::vector<std::vector<double>>().swap(gp_SATD);
	//std::vector<std::vector<double>>().swap(gp_tc);
	std::vector<std::vector<double>>().swap(gp_tc_SATD);
	// newly added
	std::vector<std::vector<double>>().swap(gp_tc_SATD_wwm);

	//***********************//
	// accuracy
	//***********************//
	// TA: calculate tau-distance on average;
	double gm_sco = .0, satd_sco = .0, rc_sco = 0.0,
		cbt_sco = 0.0, satd_wwm_sco = .0;
	for (int i = 0; i < ExperimentTimes; ++i){
		satd_sco += satd_kt_gt[i];
		gm_sco += gm_kt_gt[i];
		rc_sco += rc_kt_gt[i];
		cbt_sco += cbt_kt_gt[i];
		satd_wwm_sco += satd_wwm_kt_gt[i];
	}

	satd_sco /= (ExperimentTimes);
	gm_sco /= (ExperimentTimes);
	rc_sco /= (ExperimentTimes);
	cbt_sco /= CBT_ExperimentTimes;
	satd_wwm_sco /= (ExperimentTimes);

#if _TURE_OR_FALSE_
	of1 << "\n KanTau(SATD, GroundTruth) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << satd_kt_gt[i] << ", ";
	}
	of1 << " \nAverage KanTau(SATD, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << satd_sco << endl;


	of1 << "\n KanTau(RC, GroundTruth) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << rc_kt_gt[i] << ", ";
	}
	of1 << " \nAverage KanTau(RC, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << rc_sco << endl;


	of1 << "\n KanTau(GM, GroundTruth) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << gm_kt_gt[i] << ", ";
	}
	of1 << " \nAverage KanTau(GM, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << gm_sco << endl;

	of1 << "\n KanTau(CBT, GroundTruth) : ";
	for (int i = 0; i < CBT_ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << cbt_kt_gt[i] << ", ";
	}
	of1 << " \nAverage KanTau(CBT, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << cbt_sco << endl;

	// newly added
	of1 << "\n Truth Discovery Iteration: threshold is " << totalTime
		<< ", Stop accuracy is " << _ACCURCY << ", Execution times are:" << endl;
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << TruthIterationTimes[i] << ", ";
	}

	of1 << " \nAverage KanTau(SATD_WWM, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << satd_wwm_sco << endl;

	of1 << "\n KanTau(SATD_WWM, GroundTruth) : ";
	for (int i = 0; i < ExperimentTimes; ++i){
		of1 << std::fixed << std::setprecision(8) << satd_wwm_kt_gt[i] << ", ";
	}
	of1 << " \nAverage KanTau(SATD_WWM, GroundTruth) = "
		<< std::fixed << std::setprecision(8) << satd_wwm_sco << endl;

#endif
	SATD_tau = satd_sco;// SATD tau distance;
	RC_tau = rc_sco; // RC tau distance;
	GM_tau = gm_sco; // GM tau distance;
	CBT_tau = cbt_sco; // Crowd_BT tau distance;
	SATD_WWM_tau = satd_wwm_sco; // SATD_TT tau distance;
	

	//***********************//
	// running time
	//***********************//

	double rc_time = 0, satd_time = 0, gm_time = 0, 
		cbt_time = 0,satd_wwm_time = 0;


	// of1 << "\n SATD running time in seconds: \n";
	for (int i = 0, j = satd_elapse_time.size();
		i < j;  ++i){
		/*of1 << std::fixed << std::setprecision(8)
			<< v_satd_elapse_time[i] << ", ";*/
		satd_time += satd_elapse_time[i];
	}

	satd_time /= ExperimentTimes;

	SATD_sec = satd_time + t_voting_matrix/ExperimentTimes + SATD_WWM_TD_sec/ExperimentTimes;
		+ t_satd_tc/ExperimentTimes; // SATD running time;

	of1 << "\nSATD average time is " << std::fixed << std::setprecision(8)
		<< SATD_sec << " seconds.\n";

	of1 << "\n SATD TC step running time in seconds: \n";
	for (int i = 0; i < ExperimentTimes;  ++i){
		of1 << std::fixed << std::setprecision(8)
			<< v_tc_secs[i] << ", ";
		satd_wwm_time += satd_wwm_elapse_time[i];
	}

	of1 << "\n SATD TD step running time in seconds: \n";
	for (int i = 0;
		i < ExperimentTimes;  ++i){
		of1 << std::fixed << std::setprecision(8)
			<< v_td_secs[i] << ", ";
	}


	auto biggest = std::max_element(std::begin(v_td_secs), std::end(v_td_secs));
	double temp_max = *biggest;
	auto smallest = std::min_element(std::begin(v_td_secs), std::end(v_td_secs));
	double temp_min = *smallest;

	if (ExperimentTimes > 2)
	SATD_WWM_TD_sec = (SATD_WWM_TD_sec - temp_max - temp_min) / (ExperimentTimes -2);
	else
		SATD_WWM_TD_sec /= ExperimentTimes;
		


	biggest = std::max_element(std::begin(v_tc_secs), std::end(v_tc_secs));
	temp_max = *biggest;
	smallest = std::min_element(std::begin(v_tc_secs), std::end(v_tc_secs));
	temp_min = *smallest;

	if (ExperimentTimes > 2)
		SATD_WWM_TC_sec = (SATD_WWM_TC_sec - temp_max - temp_min) / (ExperimentTimes - 2);
	else
		SATD_WWM_TC_sec /= ExperimentTimes;
	

	satd_wwm_time /= ExperimentTimes;
	SATD_WWM_TD_sec += (t_voting_matrix / ExperimentTimes);
	SATD_WWM_Smoothing_sec /= ExperimentTimes; // smoothing time;
	SATD_WWM_HP_sec = satd_wwm_time; // find HP;

	SATD_WWM_sec = SATD_WWM_TD_sec + SATD_WWM_Smoothing_sec + SATD_WWM_TC_sec
		+ SATD_WWM_HP_sec; // SATD running time;

	of1 << "\nSATD_WWM average time is " << std::fixed << std::setprecision(8)
		<< SATD_WWM_sec << " seconds.\n";


	// of1 << "\n RC running time in seconds: \n";
	for (int i = 0, j = rc_elapse_time.size();
		i < j;  ++i){
		/*of1 << std::fixed << std::setprecision(8)
			<< v_rc_elapse_time[i] << ", ";*/
		rc_time += rc_elapse_time[i];
	}

	rc_time /= ExperimentTimes;
	RC_sec = rc_time + t_voting_matrix / ExperimentTimes; // RC running time;

	of1 << "\nRC average time is " << std::fixed << std::setprecision(8)
		<< RC_sec << " seconds.\n";

	
	// of1 << "\n GM running running time in seconds: \n";
	for (int i = 0, j = gm_elapse_time.size();
		i < j;  ++i){
		/*of1 << std::fixed << std::setprecision(8)
			<< v_gm_elapse_time[i] << ", ";*/
		gm_time += gm_elapse_time[i];
	}
	gm_time /= ExperimentTimes;
	GM_sec = gm_time + t_voting_matrix / ExperimentTimes;; // GM running time;
	of1 << "\nGM average running time is " << std::fixed << std::setprecision(8)
		<< GM_sec << " seconds.\n";


	// of1 << "\n CBT running running time in seconds: \n";
	for (int i = 0; i < CBT_ExperimentTimes;  ++i){
		/*of1 << std::fixed << std::setprecision(8)
		<< v_cbt_elapse_time[i] << ", ";*/
		cbt_time += cbt_elapse_time[i];
	}
	cbt_time /= CBT_ExperimentTimes;
	CBT_sec = cbt_time; // CBT running time;
	of1 << "\nCBT average running time is " << std::fixed << std::setprecision(8)
		<< CBT_sec << " seconds.\n";


	//release 
	vector<double>().swap(rc_kt_gt);
	vector<double>().swap(cbt_kt_gt);
	vector<double>().swap(satd_kt_gt);
	vector<double>().swap(gm_kt_gt);
	vector<double>().swap(gm_elapse_time);
	vector<double>().swap(satd_elapse_time);
	vector<double>().swap(rc_elapse_time);
	vector<double>().swap(cbt_elapse_time);
	// newly added
	//vector<double>().swap(sa_tt_elapse_time);
	vector<double>().swap(satd_wwm_elapse_time);
	//vector<double>().swap(sa_tt_kt_gt);
	vector<double>().swap(satd_wwm_kt_gt);
	//vector <int>().swap(TruthIterationTimes);
	// After you finished executing my code, **run free_all()**.
	td->free_all();
	delete td;
	return 0;
}