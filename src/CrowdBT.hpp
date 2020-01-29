/**
* @file: CrowdBT.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_CROWD_BT_H_
#define __HEADER__CROWD_SOURCING_CROWD_BT_H_
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <ctime>
#include <climits>
#include <utility> // std::pair, std::make_pair
#include <tuple>
#include <boost/math/special_functions/digamma.hpp>
#include <iostream>
#include "ktaub.hpp"
#include "MyUtility.hpp"

using namespace std;
/*
* KL(p || q) = ...
* p = N(mu_1, var_1), q = N(mu_2, var_2);
* where var_1 = sigma_1^2; and var_2 = sigma_2^2;
*/
// KL divergence of two Gaussian distributions.
// See its definition in Equation (1.113) on Page 55, PRML. M. Bishop. 
double KLGaussian(const double mu_p, const double var_p, const double mu_q,
	const double var_q);

// Beta function B(a, b) = Gamma(a)* Gamma(b)/ Gamma(a+ b);
// Note that Beat function is different from Beta distribution.
// see Wikipedia for detailed information about its definition. 
double Beta(const double & a, const double & b);

// KL divergence of two Beta distributions.
// Beta(x | a, b) = constant * x^(a -1) *(1-x)^(b -1)
// = Gamma(a + b)/ (Gamma(a)*Gamma(b)) * x^(a -1) *(1-x)^(b -1).
double KLBeta(const double alpha_p, const double beta_p,
	const double alpha_q, const double beta_q);


/* ***************************************
* for the first voting, (o_i, o_j), this pairwise comparison could generate
* the two possible voting results, i.e., o_i > o_j, and o_i < o_j.
* How will they affect the final global ranking, if they are chosen separately?
* ***************************************
*/



/* Here we want to study the problem that how two opposite outcomes of the
* first pair comparison, like (o_i, o_j) will impact the final global ranking?
* Specifically speaking,
* 1) For the first pair comparison, I do not set a specific one. Which pair
is chosen totally depends on maximum KL divergence in Equation 10 in this paper.
* 2) 1st-trial: For example, there are two cases (or two curves), one is normal,
the other is reverse. They are assigned two copies of the same initial prior distributions
(that is with the same \mu and \sigma etc, even though those parameters
are randomly generated.) Therefore, the same first pair comparison,
for example, (object 13, object 16) will be chosen for both of them.
Then they will vote (object 13 > object 16) and (object 13 < object 16),
respectively. The subsequent processes will continue as usual.

* 3) 2nd-trial: Since a different setting of initialization of those prior distributions,
probably a new pair, different from that in Fig 1, will be selected,
for example, (object 1, object 5). And so on ...
*/
void run_CrowdBT_1st_pair(const int & num_o, const int * p_num_w,
	const int num_w_size, const int & num_repeat);

// different setting of number of workers;
// meaning graphs with different degrees. 
void run_CrowdBT_diff_w_num(const int & num_o,
	const int & num_repeat);



struct INITIALS{
	double alpha_init;
	double beta_init;
	double eta_init;
	double mu_init;
	double var_init;
	double kappa_init;
};


struct PARAMS {
	// vector<double> v_alpha; // Beta distribution
	// vector<double> v_beta; // Beta distribution
	double alpha, beta, kappa; // Beta distribution
	/* the probability that the k-th annotator agrees
	 * with the true pairwise preference.
	 */
	// vector<double> v_eta;
	double eta;
	vector<double> v_mu; // mean of Gaussian;
	vector<double> v_emu;// exp(mu)
	vector<double> v_var; // variance of Gaussian;
	vector<vector<double>> vv_history;

	// 1) default constructor;
	PARAMS(){}
	// 2) user defined constructor;
	PARAMS(const int & num_o, const INITIALS & a){
		//v_alpha = vector<double>(num_o, a.alpha_init);
		//v_beta =  vector<double>(num_o, a.beta_init);
		alpha = a.alpha_init;
		beta = a.beta_init;
		kappa = a.kappa_init;
		// v_eta =   vector<double>(num_o, a.eta_init);
		eta = a.eta_init;

		srand((unsigned)time(NULL));

		v_mu = vector<double>(num_o, a.mu_init);
		v_emu = vector<double>(num_o, 0.0);


		for (int i = 0; i < num_o; ++i){
			// generate a uniform random value [0, 1]
			v_mu[i] *= rand() / double(RAND_MAX);
			v_emu[i] = exp(v_mu[i]);
		}
		v_var = vector<double>(num_o, a.var_init);
		vv_history = vector<vector<double>>(num_o, vector<double>(num_o, 0));
	}

	// 3) copy constructor;
	PARAMS(const PARAMS & p){
		//v_alpha = vector<double>(num_o, a.alpha_init);
		//v_beta =  vector<double>(num_o, a.beta_init);
		alpha = p.alpha;
		beta = p.beta;
		kappa = p.kappa;
		eta = p.eta;

		// vector : copy constructor
		v_mu = vector<double>(p.v_mu);
		v_emu = vector<double>(p.v_emu);
		v_var = vector<double>(p.v_var);
		vv_history = vector<vector<double>>(p.vv_history);
	}

	// destructor
	~PARAMS(){
		// release memory space;
		//vector<double>().swap(v_alpha);
		//vector<double>().swap(v_beta);
		//vector<double>().swap(v_eta);
		vector<double>().swap(v_mu);
		vector<double>().swap(v_emu);
		vector<double>().swap(v_var);
		vector<vector<double>>().swap(vv_history);
	}

};

struct NEW_PARAMS {
	// we assume object i > object j;
	// for the opposite case (i.e., object j > object i), 
	// just swap i and j before calling this function,
	// or consider always "mu_i_new" as the preferred object,
	// and "mu_j_new" as the less one.
	double mu_i_new,
		mu_j_new,
		var_i_new,
		var_j_new,
		alpha_new,
		beta_new;
};




class Crowd_BT{

private:
	int num_o; // number of objects;
	PARAMS params;
public:
	// constructor
	Crowd_BT(const int & num_o, const INITIALS & a){
		params = PARAMS(num_o, a);
		this->num_o = num_o;
	}

	// copy constructor
	// constructor
	Crowd_BT(const Crowd_BT & bt){
		params = PARAMS(bt.params);
		this->num_o = bt.num_o;
	}
	// destructor
	~Crowd_BT(){ params.~PARAMS(); }

	// To rank objects by sorting the obtained {mu_i}(means of Gaussians.);
	// return the object ranking in an ascending order.
	int * get_ascending_ranking();
	/*Equation (19), (16)*/
	void get_C(const int & i, const int & j, double & c1, double & c);

	
	/*
	* Update the parameters based on the equations shown in the paper:
	* "Pairwise Ranking Aggregation in a CrowdSourced Setting".
	* The equation will be named according to the number originally shown
	* in this paper.
	* Pay attention to this function:
	* its first parameter means the preferred object,
	* therefore, if o_j > o_i, we should input j as the first
	* parameter to this function. For the opposite case of o_i > o_j, 
	* we should pass i as the first parameter to this function.
	*/
	NEW_PARAMS get_updated_parameters(const int & i, const int& j,
		double & c1, double & c);


	/* v_layers describes how to measure the similarity of objects to be ranked.
	 * If we have L = 20 objects, the ranking is 1, 2, ..., L.
	 * We can set the ambiguity level, for example, we set v_layers = [0.25, 0.5, 1],
	 * generating the following ranking layers:
	 *  1) level 1: [0, 0.25]* L  = [0, 5];
	 *  2) level 2: [0.25, 0.75] * L = [5, 14];
	 *  3) level 3: [0.75, 1.0] * L = [14, 20];
	 * If object 1 and object 2 come from the same layer, like, they belong to level 2,
	 * there is a possibility that wrong ranking will occur,
	 * according to some probability distribution.
	 * We will implement the above idea in the following function.
	 */
	bool get_preference(const vector<double> & v_scores,
		const int & i, const int & j, const vector<double> & v_layers);

	/* Select a pair (o_i, o_j) for some annotator k which
	 * maximize the expected information gain in E.q.(10).
	 * parameter: gamma, exploration-exploitation trade-off,
	 * see Section 5.1.2 in this paper.
	 */
	//**********************************
	// before actual voting, try to select a pair of objects o_i and o_j,
	// such that they can maximize the expected information gain 
	// in E.q.(10) in this paper.
	//**********************************
	std::tuple<int, int, double> get_pair(const double & gamma);

	/* assign the new parameters to the member variables in this class.
	 * Inputs:
	 * - i_max: (object i, object j) with maximum KL divergence;
	 * - j_max: (object i, object j) with maximum KL divergence;
	 */
	void keep_new_params(const NEW_PARAMS & new_para, const int & i_max,
		const int & j_max);


	/* For one worker, like worker k, he/she will select a pair (o_i, o_j),
	* 1) after he/she do once (i.e., 1 time ) voting, we could update the 
	* parameters of those distribution according to Bayes' Theorem.
	* 2) If we consider the wrong voting probability that a worker will give wrong
	* voting result, probably due to the high similarity between o_i, o_j, or
	* due to the work quality or even the malicious behavior of this worker.
	*/
	vector<pair<int, int>> repeated_voting_1_worker(const int * p_rank_gt,
		const int & num_repeat, const vector<double> & v_scores,
		const vector<double> & v_layers, const double & gamma,
		const bool & IsDisplay, const bool & IsReverse,
		const pair<int, int> & p_ij);


	// calculate the Kendall-tau distance between two ranking:
	//  1) ranking 1: is the input, usually the ground truth ranking or
	//     the result obtained from other algorithms;
	//  2) ranking 2: is the ranking obtained by this algorithm. 
	double get_kendall_tau_distance(
		// ranking 1
		const int * p_rank_gt, /*ground truth ranking*/
		const bool & IsDisplay);

	void Crowd_BT::repeated_voting_1_worker(
		const vector<int> & v_scores, // scores of each object.
		const int & n, // vertex number
		const int & distribution,
		const int & num_repeat,
		//const double & stddev, // Gaussian variance to control the worker's quality.
		const int& quality, // different levels of error-rate for the worker's quality.
		const bool & IsDisplay, const pair<int, int> & p_ij,
		const bool & Is_Beta_Simu_data_type // == true, using Beta distribution to simulate the voting result.
		);

};

double one_worker_reverse_or_not(const int & num_o, const int num_w,
	const int * p_rank_gt,
	const int & num_repeat, const vector<double> & v_scores,
	const vector<double> & v_layers,
	const double & gamma, const bool & IsDisplay, const bool & IsReverse,
	Crowd_BT & c_bt1);

#endif