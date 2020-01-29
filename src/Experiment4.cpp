/**
* @file: Experiment4.cpp
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

extern double beta_a;
extern double beta_b;

void run_CrowdBT_for_SA_Comparison(
	double & tau,
	double & sec,
	const int & vertexNum,
	// pay attention here:
	// the order is i --> j, means to prefer i than j, is different from SA;.
	// make sure ascending order!!
	const int * AscendingGTr, // ground truth ranking, ascending order based on their scores.
	const int & distribuion_type,
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & w, // number of workers;
	const int & num_repeat,
	const bool & Is_Beta_Simu_data_type // == true, using Beta distribution to simulate the voting result.
	){
	
	
	std::srand(unsigned(std::time(0)));
	vector<pair<int, int>> v_voting;
	auto t2 = std::chrono::high_resolution_clock::now();
	//******************
	// initial values;
	//******************
	INITIALS a;
	a.alpha_init = beta_a; // e.g. , == 10;
	a.beta_init = beta_b; //e.g. ==  5;
	a.eta_init = 1;
	a.mu_init = 1;
	a.var_init = 1.0 / 9.0;
	a.kappa_init = 1.0e-4;
	// make sure the ascending order of the elements.
	vector<int> v_scores(vertexNum, 0);
	// 5 <= gamma <= 10 for experiments;
	double gamma = 5;
	// ranking scores, the larger the better.
	// the scores are generated based on the ranking ground truth.
	for (int i = 0; i < vertexNum; ++i)
		v_scores[AscendingGTr[i]] = i+1;

	
	const bool IsDisplay = false;
	Crowd_BT c_bt1(vertexNum, a);

	for (int ww = 0; ww < w; ++ww){ // each worker

	tuple<int, int, double> t_i_j = c_bt1.get_pair(gamma);
	int i = std::get<0>(t_i_j), j = std::get<1>(t_i_j);
	double max_kl_div = std::get<2>(t_i_j);

	
	// Let him/her repeatly vote for the pair-wise comparison;
	// i..e, (o_i, o_j) , due to the wrong answer probability, especially 
	// for very similar objects.
	c_bt1.repeated_voting_1_worker(
		v_scores, vertexNum, distribuion_type, num_repeat, quality,
		IsDisplay, make_pair(i, j), Is_Beta_Simu_data_type);

	}
	// final ranking after all the workers' voting.
    // get ascending order based on their scores, that is the means of Gaussian distributions.
	tau = c_bt1.get_kendall_tau_distance(AscendingGTr, IsDisplay);
	auto t3 = std::chrono::high_resolution_clock::now();
	auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
	sec = dt;
}