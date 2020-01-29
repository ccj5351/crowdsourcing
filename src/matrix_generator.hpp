/**
* @file: matrix_generator.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

// some of the codes are generated by Boxiang Do.
#ifndef matrix_generator_hpp
#define matrix_generator_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include "truth_discovery.hpp"
#include <random>
#include <boost/math/special_functions/beta.hpp>

#define UNIFORM_DISTRIBUTION 0
#define NORMAL_DISTRIBUTION 1
#define NEARLY_ZERO 1.0E-5

#define LOW_QUALITY 0
#define MEDIUM_QUALITY 1
#define HIGH_QUALITY 2

#define low_guassian 1
#define medium_guassian  0.1
#define high_guassian  0.01
#define high_uniform_left  0
#define high_uniform_right 0.2
#define medium_uniform_left  0.1
#define medium_uniform_right 0.3
#define low_uniform_left  0.2
#define low_uniform_right  0.4

/*
double low_guassian = 1;
double medium_guassian = 0.1;
double high_guassian = 0.01;
double high_uniform_left = 0;
double high_uniform_right = 0.2;
double medium_uniform_left = 0.1;
double medium_uniform_right = 0.3;
double low_uniform_left = 0.2;
double low_uniform_right = 0.4;
*/


using namespace std;
bool generate_sim_vote_for_one_task_for_CrowdBT(
	const int & n, // vertex number
	const int & rank_o_i,
	const int & rank_o_j, // less rank is better.
	const double & a, //Reta distribution parameter.
	const double & b //Beta distribution parameter.
	);

bool generate_sim_vote_for_one_task(
	const int & n, // vertex number
	const int & rank_o_i,
	const int & rank_o_j, // less rank is better.
	const double& stddev /*Gaussian variance to control the worker's quality.*/
	);


vector<vector<int>> generate_sim_vote(const int& n, bool **GT,
	// ground_truth: <key = object_i, value = rank of object_i>;
	const map<int, int>& ground_truth,
	const double& stddev, /*Gaussian variance to control the worker's quality.*/
	const int & max_per_worker // maximum number of pairwise comparison per worker can do.
	// if max_per_worker < 0, that means doing all the pairwise comparasions.
	);




// for all the workers, gathering their voting results.
vector<vector<vector<int>>> generate_sim_votes(
	const int& n, const int& w, bool **GT,
	const map<int, int>& ground_truth, const int& distribution,
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & max_per_worker /*maximum number of pairwise comparison per worker can do.*/
	);



void generate_sim_votes_taskID_who(
	const int& n, const int& w, bool **GT,
	const vector<vector<vector<int>>> & matrice,
	// store the tasks;
	// The task-IDs' are assigned based on way the double-for-loop is executed.
	vector<vector <td_wieghts::client_vote>> & v_client_votes
	);


vector<vector <td_wieghts::client_vote>>
generate_aggregated_matrix_SATD(bool **GT,
	const vector<vector<vector<int>>> & matrice,
	std::vector<std::vector<double>> & gp, 
	vector<double> & worker_weights,
	td_wieghts::truth_discovery * td,
    const int & max_task_num,
	int & count, // Truth discovery execution times.
	/*Smothing_Type : */
	const int & smothingType = DIRECT_USE_WEIGHTS
	);

double small_prob(const double & worker_weight, const double & mean,
	const double & prob_threshold, /*a threshold of the small_prob to be generated.*/
	const double & stddev_scalar);

void Smoothing_UpdateProb_with_Workers_Quality(
	const vector<double> & worker_weights,
	bool **GT, const int & n, /*object number*/
	const double & prob_threshold, /*a threshold of the small_prob to be generated.*/
	const double & mean, const double & stddev_scalar,
	const vector<vector <td_wieghts::client_vote>> & v_client_votes,
	vector<vector<double>> & gtc);

/* You have to specify the default values
* for the arguments only in the declaration
* but not in the definition.*/
std::vector<std::vector<double> > generate_aggregated_matrix(
	const int& n, const int& w, bool **GT,
	const std::map<int, int>& ground_truth, const int& distribution,
	const int& quality, // different levels of error-rate for the worker's quality.
	const int & max_per_worker, // maximum number of pairwise comparison per worker can do.
	/*Smothing_Type : */
	const int & Smothing_Type = DIRECT_USE_WEIGHTS
	);


std::vector<std::vector<double> >
generate_aggregated_matrix(
const std::vector<std::vector<std::vector<int>>> & matrice,
/*Smothing_Type : */
const int & smothingType = DIRECT_USE_WEIGHTS
);

// added by CCJ
std::vector<std::vector<int> >
generate_aggregated_voting_matrix(const vector<vector<vector<int>>> &  matrice);



#endif /* matrix_generator_hpp */