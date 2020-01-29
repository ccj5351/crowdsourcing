/*
* @file: task_assignment.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "CrowdBT.hpp"
#include "matrix_generator.hpp"

/* Here I implemented the algorithm in the paper:  * Pairwise Ranking Aggregation in a Crowdsourced Setting. */

using namespace std;
extern std::default_random_engine generator_main;
/*
* KL(p || q) = ...
* p = N(mu_1, var_1), q = N(mu_2, var_2);
* where var_1 = sigma_1^2; and var_2 = sigma_2^2;
*/
// KL divergence of two Gaussian distributions.
// See its definition in Equation (1.113) on Page 55, PRML. M. Bishop. 
double KLGaussian(const double mu_p, const double var_p, const double mu_q,
	const double var_q){
	// log Returns the natural logarithm of x.
	double kl = log(var_q / var_p) + var_p / var_q +
		(mu_p - mu_q)*(mu_p - mu_q) / var_q - 1;
	return kl / 2;
}

// Beta function B(a, b) = Gamma(a)* Gamma(b)/ Gamma(a+ b);
// Note that Beat function is different from Beta distribution.
// see Wikipedia for detailed information about its definition. 
double Beta(const double & a, const double & b){
	return tgamma(a) * tgamma(b) / tgamma(a + b);
}

// KL divergence of two Beta distributions.
// Beta(x | a, b) = constant * x^(a -1) *(1-x)^(b -1)
// = Gamma(a + b)/ (Gamma(a)*Gamma(b)) * x^(a -1) *(1-x)^(b -1).
double KLBeta(const double alpha_p, const double beta_p,
	const double alpha_q, const double beta_q){
	// Returns the gamma function of x.
	double B_p = Beta(alpha_p, beta_p);
	double B_q = Beta(alpha_q, beta_q);
	double phi_a_p = boost::math::digamma(alpha_p);
	double phi_b_p = boost::math::digamma(beta_p);
	double phi_ab_p = boost::math::digamma(alpha_p + beta_p);
	double kl = log(B_q / B_p) + (alpha_p - alpha_q) * phi_a_p +
		(beta_p - beta_q) * phi_b_p + (alpha_q - alpha_p + beta_q - beta_p)*phi_ab_p;
	return kl;
}

/* ***************************************
* for the first voting, (o_i, o_j), this pairwise comparison could generate
* the two possible voting results, i.e., o_i > o_j, and o_i < o_j.
* How will they affect the final global ranking, if they are chosen separately?
* ***************************************
*/
double one_worker_reverse_or_not(const int & num_o, const int num_w,
	const int * p_rank_gt,
	const int & num_repeat, const vector<double> & v_scores,
	const vector<double> & v_layers,
	const double & gamma, const bool & IsDisplay, const bool & IsReverse,
	Crowd_BT & c_bt1){


	//**********************************
	// For worker k = 0, i.e., the first worker.
	// before his/her voting, try to select a pair of objects o_i and o_j,
	// such that they can maximize the expected information gain 
	// in E.q.(10) in this paper.
	//**********************************

	tuple<int, int, double> t_i_j = c_bt1.get_pair(gamma);
	int i = std::get<0>(t_i_j), j = std::get<1>(t_i_j);
	double max_kl_div = std::get<2>(t_i_j);


	//**********************************
	// For the first voting, (o_i, o_j), this pairwise comparison could generate 
	// the two possible voting results, i.e., o_i > o_j, and o_i < o_j.
	//**********************************
	int k = 0;
	if (IsDisplay){
		std::cout << "Worker " << k << " chooses pairwise comparison ";
		std::cout << "( " << i << ", " << j
			<< " ), due to its maximum KL divergence = "
			<< max_kl_div << std::endl;
	}

	// Let him/her repeatly vote for the pair-wise comparison;
	// i..e, (o_i, o_j) , due to the wrong answer probability, especially 
	// for very similar objects.
	// Pay attention to this parameter "IsReverse" when calling the following function;
	vector<pair<int, int>> v_votes = c_bt1.repeated_voting_1_worker(p_rank_gt, num_repeat,
		v_scores, v_layers, gamma, IsDisplay, IsReverse, make_pair(i, j));


	for (int k = 1; k < num_w; ++k){/* for worker k.*/

		//**********************************
		// before voting, try to select a pair of objects o_i and o_j,
		// such that they can maximize the expected information gain 
		// in E.q.(10) in this paper.
		//**********************************

		tuple<int, int, double> t_i_j = c_bt1.get_pair(gamma);
		int i = std::get<0>(t_i_j), j = std::get<1>(t_i_j);
		double max_kl_div = std::get<2>(t_i_j);
		if (IsDisplay){
			std::cout << "Worker " << k << " chooses pairwise comparison ";
			std::cout << "( " << i << ", " << j
				<< " ), due to its maximum KL divergence = "
				<< max_kl_div << std::endl;
		}

		// for annotator k, let him/her repeatly vote for the pair-wise comparison;
		// i..e, (o_i, o_j) , due to the wrong answer probability, especially 
		// for similar objects.
		// Pay attention to this parameter "IsReverse = false" here,
		// when calling the following function.
		vector<pair<int, int>> v_votes = c_bt1.repeated_voting_1_worker(
			p_rank_gt, num_repeat, v_scores, v_layers, gamma, IsDisplay,
			false, make_pair(i, j));

		// c_bt1.get_kendall_tau_distance(p_rank_gt, IsDisplay);

	} /*end of worker k.*/

	// final ranking after all the workers' voting.
	double tau = c_bt1.get_kendall_tau_distance(p_rank_gt, IsDisplay);
	return tau;
}


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
	const int num_w_size, const int & num_repeat){

	std::srand(unsigned(std::time(0)));
	vector<pair<int, int>> v_voting;
	//******************
	// initial values;
	//******************
	INITIALS a;
	a.alpha_init = 10;
	a.beta_init = 1;
	a.eta_init = 1;
	a.mu_init = 1;
	a.var_init = 1.0 / 9.0;
	a.kappa_init = 1.0e-4;
	// make sure the ascending order of the elements.
	vector<double> v_layers = { 0.25, 0.75, 1.0 };
	vector<double> v_scores;
	// 5 <= gamma <= 10 for experiments;
	double gamma = 5;
	// ranking scores, randomly generated.
	for (int i = 0; i < num_o; ++i)
		v_scores.push_back(i + 1);
	// using built-in random generator:
	std::random_shuffle(v_scores.begin(), v_scores.end());

	// <object i's index, object i 's score>
	vector<std::pair<int, double>> v_pair_rank_gt(num_o);
	vector<std::pair<int, double>> v_pair_rank_crowd_bt(num_o);
	for (int i = 0; i < num_o; ++i){
		v_pair_rank_gt[i] = make_pair(i, v_scores[i]);
	}
	// sorting based on the scores.
	std::sort(v_pair_rank_gt.begin(), v_pair_rank_gt.end(),
		[](const std::pair<int, double> & arg1,
		const std::pair<int, double> & arg2){return arg1.second <= arg2.second; });

	int * p_rank_gt = new int[num_o];
	for (int i = 0; i < num_o; ++i){
		p_rank_gt[i] = v_pair_rank_gt[i].first;
	}

	const bool IsDisplay = false;
	Crowd_BT c_bt0(num_o, a);
	// kendall tau distances;
	vector<pair<double, double>> v_p_taus;
	// for different degrees of the graph;
	for (int ww = 0; ww < num_w_size; ++ww){
		const int num_w = p_num_w[ww];
		Crowd_BT c_bt1(c_bt0), c_bt2(c_bt0);// two copies;
		//****************
		// Case 1: worker 1, 
		// normal first-task voting;
		//****************
		bool IsReverse = false;
		double tau1 = one_worker_reverse_or_not(num_o, num_w, p_rank_gt, num_repeat,
			v_scores, v_layers, gamma, IsDisplay, IsReverse, c_bt1);


		//****************
		// Case 2: worker 1, 
		// reverse first-task voting;
		//****************
		IsReverse = true;
		double tau2 = one_worker_reverse_or_not(num_o, num_w, p_rank_gt, num_repeat,
			v_scores, v_layers, gamma, IsDisplay, IsReverse, c_bt2);

		cout << "Case 1 (Normal) : final Kendall Tau = " << tau1 << endl
			<< "Case 2 (Reverse): final Kendall Tau = " << tau2 << endl;
		v_p_taus.push_back(make_pair(tau1, tau2));
	}/*end of different degrees of the graph.*/

	if (IsDisplay){
		for (int i = 0; i < v_p_taus.size(); ++i){
			std::cout << " degree = " << p_num_w[i] << ", tau distance = "
				<< v_p_taus[i].first << " (normal) , " << v_p_taus[i].second
				<< " (reverse).\n";
		}
	}
	/*
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	string cur_time = to_string((now->tm_mon + 1)) + "-"
	+ to_string(now->tm_mday) + "-" + to_string(now->tm_hour)
	+ "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec);
	*/

	string s_ofile = "E:/OpenCVProjects_CCJ/CrowdSourcing2/Experiments/Crowd_BT/object-"
		+ to_string(num_o) + ".txt";
	std::ofstream  of23;
	of23.open(s_ofile.c_str(), std::ofstream::out);
	if (!of23.is_open())
		std::cout << "ofstream cannot open " << s_ofile << endl;

	for (int i = 0; i < v_p_taus.size(); ++i){
		of23 << " degree = " << p_num_w[i] << ", tau distance = "
			<< v_p_taus[i].first << " (normal),  " << v_p_taus[i].second
			<< " (reverse)." << endl;
	}

	/*
	of23 << endl << endl;
	for (int i = 0; i < v_p_taus.size(); ++i){
	of23 << p_num_w[i] << "\t" << v_p_taus[i].first << "\t" << v_p_taus[i].second << endl;
	}*/
	// outputing Matlab code;
	of23 << " x = 2:1: " << num_o - 1 << ";" << endl << " y1 = [ ";
	for (int i = 0; i < v_p_taus.size(); ++i){
		of23 << v_p_taus[i].first << " , ";
	}
	of23 << "];\n";
	of23 << "y2 = [ ";
	for (int i = 0; i < v_p_taus.size(); ++i){
		of23 << v_p_taus[i].second << " , ";
	}
	of23 << "];\n f1 = figure;\n plot(x, y1, x, y2);\n title(' object num = "
		<< num_o << "');\n"
		<< " xlabel(' numbers of worker');\n ylabel('Kendall-Tau distance');\n "
		<< "legend('normal', 'reverse');\n saveas(f1, 'object-" << num_o
		<< ".jpg');" << endl;
	of23.close();
}



// different setting of number of workers;
// meaning graphs with different degrees. 
void run_CrowdBT_diff_w_num(const int & num_o,
	const int & num_repeat){
	int num_w_size = num_o - 2;
	int * p_num_w = new int[num_w_size];
	// different number of workers, 
	// that is to say, different workers means different degrees in the
	// corresponding task assignment graphs.
	for (int i = 0; i < num_w_size; ++i){
		p_num_w[i] = std::ceil((double)(num_o * (2 + i)) /
			(double)(2 * num_repeat));
	}
	run_CrowdBT_1st_pair(num_o, p_num_w, num_w_size, num_repeat);
	delete[] p_num_w;
}


// To rank objects by sorting the obtained {mu_i}(means of Gaussians.);
// return the object ranking in an ascending order.
int * Crowd_BT::get_ascending_ranking(){
		int * p_rank_crowd_bt = new int[num_o];
		// <object i's index, object i 's score>
		vector<std::pair<int, double>> v_pair_rank_crowd_bt(num_o);
		for (int i = 0; i < num_o; ++i){
			v_pair_rank_crowd_bt[i] = make_pair(i, params.v_mu[i]);
		}

		std::sort(v_pair_rank_crowd_bt.begin(), v_pair_rank_crowd_bt.end(),
			[](const std::pair<int, double> & arg1, 
			const std::pair<int, double> & arg2)
		{return arg1.second <= arg2.second; });


		for (int i = 0; i < num_o; ++i){
			p_rank_crowd_bt[i] = v_pair_rank_crowd_bt[i].first;
		}

		// release memory;
		vector<std::pair<int, double>>().swap(v_pair_rank_crowd_bt);
		return p_rank_crowd_bt;
	}

/*Equation (19), (16)*/
void Crowd_BT::get_C(const int & i, const int & j, double & c1, double & c){
		// Eq (19), to get C1;
		// exp(mu_i) = params.v_emu[i], 
		// exp(mu_j) = params.v_emu[j].

		double denomi = pow((params.v_emu[j] + params.v_emu[i]), 3);
		// c1
		c1 = params.v_emu[i] / (params.v_emu[i] + params.v_emu[j])
			+ 0.5*(params.v_var[i] + params.v_var[j])*params.v_emu[i] *
			params.v_emu[j] * (params.v_emu[j] - params.v_emu[i]) / denomi;
		// c
		c = (c1 * params.alpha + (1 - c1) * params.beta)
			/ (params.alpha + params.beta);
	}

	
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
NEW_PARAMS Crowd_BT::get_updated_parameters(const int & i, const int& j,
		double & c1, double & c){

		NEW_PARAMS new_para;

		// Eq (12) & Eq (13);
		double term1 = params.alpha *params.v_emu[i] / (params.alpha *
			params.v_emu[i] + params.beta* params.v_emu[j]) - params.v_emu[i]
			/ (params.v_emu[i] + params.v_emu[j]);
		new_para.mu_i_new = params.v_mu[i] + params.v_var[i] * term1;// Eq.(12);
		new_para.mu_j_new = params.v_mu[j] - params.v_var[j] * term1;// Eq.(13);

		double term2 = (params.alpha * params.v_emu[i] * params.beta* params.v_emu[j])
			/ pow((params.alpha * params.v_emu[i] + params.beta* params.v_emu[j]), 2)
			- ((params.v_emu[i] * params.v_emu[j])) 
			/ pow((params.v_emu[i] + params.v_emu[j]), 2);

		// Eq. (14);
		double temp = 1 + params.v_var[i] * term2;
		new_para.var_i_new = params.v_var[i] * std::max(temp, params.kappa);
		// Eq. (15);
		temp = 1 + params.v_var[j] * term2;
		new_para.var_j_new = params.v_var[j] * std::max(temp, params.kappa);

		//Preparing for calculating E.q. (17) and (18):
		get_C(i, j, c1, c);
		// expectation of eta;
		double E_eta = (c1*(params.alpha + 1)*params.alpha +
			(1 - c1)* params.alpha*params.beta) / (c*(params.alpha + params.beta + 1)*
			(params.alpha + params.beta));
		// expectation of eta^2;
		double E_eta2 = (c1*(params.alpha + 2)*(params.alpha + 1)*params.alpha +
			(1 - c1)* (params.alpha + 1)*params.alpha*params.beta) / (c*(params.alpha + params.beta + 2)
			*(params.alpha + params.beta + 1)* (params.alpha + params.beta));
		//E.q. (17)
		new_para.alpha_new = (E_eta - E_eta2)*E_eta / (E_eta2 - E_eta * E_eta);
		//E.q. (18)
		new_para.beta_new = (E_eta - E_eta2)* (1 - E_eta) / (E_eta2 - E_eta * E_eta);
		return new_para;
	}


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
bool Crowd_BT::get_preference(const vector<double> & v_scores,
		const int & i, const int & j, const vector<double> & v_layers){
		double score_i = v_scores[i], score_j = v_scores[j];
		int level_i = 0, level_j = 0;
		for (int k = 0; k < v_layers.size(); ++k){
			// make sure: the input parameter "v_layers" is sorted in the 
			// ascending order before.
			// Also make sure the range of score values (exactly speaking, 
			// it should be the ranking index), is in [0, L],
			// where L is the number of the objects for ranking.
			if (score_i > v_layers[k])
				level_i = k;
			if (score_j > v_layers[k])
				level_j = k;
		}

		bool Is_i_Better = (score_i > score_j);
		// If object i and object j come from the same layer,
		// there is a possibility that wrong ranking result will occur.
		if (level_i == level_j){
			double p = exp(-(score_i + score_j)); // probability of wrong answer.
			srand((unsigned)time(NULL));
			double r = rand() / double(RAND_MAX);
			if (r < p)
				Is_i_Better = !(Is_i_Better);
		}
		return Is_i_Better;
	}

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
std::tuple<int, int, double> Crowd_BT::get_pair(const double & gamma){
		// Query the annotator k on the preference between o_i and o_j.
		NEW_PARAMS new_para;
		double c1 = .0, c = .0,
			max_kl_div = -1e9; // - 1 * 10^9;
		int max_kl_div_idx_i = 0, max_kl_div_idx_j = 0;
		for (int i = 0; i < num_o; ++i)
			for (int j = 0; j < num_o; ++j)
			{
				if (i != j){
					// case 1: assuming o_i > o_j, 
					// give this new observation - o_i > o_j, 
					// then to update the parameters in the posterior.
					new_para = get_updated_parameters(i, j, c1, c);
					// KL(p||q) = KL(p = N(mu_new, var_new) || 
					//               q = N(mu_old, var_old));

					// KL_oi_ij: means KL divergence between the posterior 
					// and prior for object i, given the observation o_i > o_j,
					// that is, the annotator k prefers o_i than o_j. 
					double KL_oi_ij = KLGaussian(new_para.mu_i_new,
						new_para.var_i_new, params.v_mu[i], params.v_var[i]);

					// KL_oj_ij: means KL divergence between the posterior 
					// and prior for object j, given the observation o_i > o_j,
					// that is, the annotator k prefers o_i than o_j.
					double KL_oj_ij = KLGaussian(new_para.mu_j_new,
						new_para.var_j_new, params.v_mu[j], params.v_var[j]);

					// KL_wk_ij: means KL divergence between the posterior 
					// and prior for annotator k, given the observation o_i > o_j,
					// that is, the annotator k prefers o_i than o_j.
					double KL_wk_ij = KLBeta(new_para.alpha_new, new_para.beta_new, 
						params.alpha, params.beta);

					double KL_div_1 = c*(KL_oi_ij + KL_oj_ij + gamma * KL_wk_ij);

					// case 2: assuming assuming o_j > o_i
					// give this new observation - o_j > o_i, then to update 
					// the parameters in the posterior.
					// pay attention to the following function get_updated_parameters:
					// its first parameter means the preferred object,
					// therefore, if o_j > o_i, we should input j as the first
					// parameter to this function.
					new_para = get_updated_parameters(j, i, c1, c);
					// swap, since mu_i should decrease, but mu_j should increase.
					swap(new_para.mu_i_new, new_para.mu_j_new);
					swap(new_para.var_i_new, new_para.var_j_new);
					// KL_oj_ji: means KL divergence between the posterior 
					// and prior for object j, given the observation o_j > o_i,
					// that is, the annotator k prefers o_j than o_i. 
					double KL_oi_ji = KLGaussian(new_para.mu_i_new,
						new_para.var_i_new, params.v_mu[i], params.v_var[i]);
					double KL_oj_ji = KLGaussian(new_para.mu_j_new,
						new_para.var_j_new, params.v_mu[j], params.v_var[j]);
					double KL_wk_ji = KLBeta(new_para.alpha_new, new_para.beta_new, 
						params.alpha, params.beta);
					// E.q. (10)
					double KL_div_2 = c*(KL_oi_ji + KL_oj_ji + gamma * KL_wk_ji);
					double KL_div = KL_div_1 + KL_div_2;
					if (max_kl_div < KL_div){
						// keep the currently large KL divergence, and the indices i and j.
						max_kl_div = KL_div;
						max_kl_div_idx_i = i;
						max_kl_div_idx_j = j;
					}

				} /*end of object j*/
			}/*end of object i*/

		return std::make_tuple(max_kl_div_idx_i, max_kl_div_idx_j, max_kl_div);
	}

	/* assign the new parameters to the member variables in this class.
	 * Inputs:
	 * - i_max: (object i, object j) with maximum KL divergence;
	 * - j_max: (object i, object j) with maximum KL divergence;
	 */
void Crowd_BT::keep_new_params(const NEW_PARAMS & new_para, const int & i_max,
		const int & j_max){
		// object i;
		params.v_mu[i_max] = new_para.mu_i_new;
		params.v_var[i_max] = new_para.var_i_new;
		params.v_emu[i_max] = exp(new_para.mu_i_new);
		// object j;
		params.v_mu[j_max] = new_para.mu_j_new;
		params.v_var[j_max] = new_para.var_j_new;
		// annotator k;
		params.alpha = new_para.alpha_new;
		params.beta = new_para.beta_new;
	}


	/* For one worker, like worker k, he/she will select a pair (o_i, o_j),
	* 1) after he/she do once (i.e., 1 time ) voting, we could update the 
	* parameters of those distribution according to Bayes' Theorem.
	* 2) If we consider the wrong voting probability that a worker will give wrong
	* voting result, probably due to the high similarity between o_i, o_j, or
	* due to the work quality or even the malicious behavior of this worker.
	*/
vector<pair<int, int>> Crowd_BT::repeated_voting_1_worker(
	const int * p_rank_gt,
		const int & num_repeat, const vector<double> & v_scores, 
		const vector<double> & v_layers, const double & gamma, 
		const bool & IsDisplay, const bool & IsReverse, 
		const pair<int, int> & p_ij){


		// Repeated voting begins here;
		// As for how many times he/she can repeatedly vote, depends
		// totally on the budget predefined.
		vector<pair<int, int>> v_votes;
		int i = p_ij.first, j = p_ij.second;
		if (IsDisplay){
			std::cout << "  Voting begins:\n";
		}
		// for annotator k, let him repeatly vote for the pair-wise comparison;
		// i.e., (o_i, o_j) , due to the wrong answer probability, especially 
		// for similar objects.
		for (int r = 0; r < num_repeat; ++r){
			bool pref = this->get_preference(v_scores, i, j, v_layers);
			// the "IsReverse" parameter is used to compare the impacts to 
			// the final global ranking,
			// incurred by two opposites voting results, that is:
			//    1) (o_i > o_j); 
			//    2) (o_i < o_j).
			if (IsReverse){
				if (!pref){ //  (o_i < o_j);
					if (IsDisplay){
						// due to IsReverse = True in this case.
						cout << " ----  " << r << " voting : " << " o_" << i
							<< " > o_" << j << ", ";
					}
					// make sure "i" means the better one;
					std::swap(i, j);
				}
				else{
					if (IsDisplay){
						cout << " ----  " << r << " voting : " << " o_" << j
							<< " > o_" << i << ", ";
					}
				}
				std::swap(i, j);
			}
			
			else{
				if (!pref){ //  (o_i < o_j);
					if (IsDisplay){
						cout << " ----  " << r << " voting : " << " o_" << j
							<< " > o_" << i << ", ";
					}
					// make sure "i" means the better one;
					std::swap(i, j);
				}
				else{
					if (IsDisplay){
						cout << " ----  " << r << " voting : " << " o_" << i
							<< " > o_" << j << ", ";
					}
				}
			}

			
			// store the voting results by worker k;
			v_votes.push_back(make_pair(i, j));

			//**********************************
			// once the worker finishes voting to generate a preference result,
			// for example (o_j > o_i). Then based on this new observation,
			// we update the parameters according to
			// E.q. (12),(13),(14),(15),(17) and (18).
			//**********************************
			double c1, c;
			NEW_PARAMS new_para = this->get_updated_parameters(i, j, c1, c);

			// update parameters 
			this->keep_new_params(new_para, i, j);


			// new ranking
			int * p_rank_crowd_bt = this->get_ascending_ranking();
			// Kendall-tau distance between two rankings;
			double tau = kendallSmallN(p_rank_gt, "GT", p_rank_crowd_bt, 
				"Crowd_BT", num_o, false);
			if (IsDisplay){
				std::cout << "Kendall-Tau distance is " << tau << endl;
			}
			// release memory;
			delete[] p_rank_crowd_bt;
		} /*end of repeated voting for one worker.*/

		return v_votes;
	}


// As baseline for experiment comparison.
void Crowd_BT::repeated_voting_1_worker(
	const vector<int> & v_scores, // scores of each object.
	const int & n, // vertex number
	const int & distribution,
	const int & num_repeat, 
	//const int & mean, // Gaussian mean to control the worker's quality.
	// const double & stddev, // Gaussian variance to control the worker's quality.
	const int& quality, // different levels of error-rate for the worker's quality.
	const bool & IsDisplay, const pair<int, int> & p_ij,
	const bool & Is_Beta_Simu_data_type // == true, using Beta distribution to simulate the voting result.
	){


	// Repeated voting begins here;
	// As for how many times he/she can repeatedly vote, depends
	// totally on the budget predefined.
	int i = p_ij.first, j = p_ij.second;
	// make sure the ascending order of the elements.
	// ranking scores, randomly generated.
	int score_i = v_scores[i], score_j = v_scores[j];
	int rank_i = n - score_i, rank_j = n - score_j;
	

	double stddev = 0.001;
	std::normal_distribution<double> normal_dist;
	std::uniform_real_distribution<double> uniform_dist;

	if (distribution == NORMAL_DISTRIBUTION) {
		switch (quality) {
		case LOW_QUALITY:
			normal_dist = std::normal_distribution<double>(0, low_guassian);
			break;
		case MEDIUM_QUALITY:
			normal_dist = std::normal_distribution<double>(0, medium_guassian);
			break;
		case HIGH_QUALITY:
			normal_dist = std::normal_distribution<double>(0, high_guassian);
			break;
		default:
			break;
		}

		stddev = abs(normal_dist(generator_main));
	}
	else {
		switch (quality) {
		case LOW_QUALITY:
			uniform_dist = std::uniform_real_distribution<double>(low_uniform_left, low_uniform_right);
			break;
		case MEDIUM_QUALITY:
			uniform_dist = std::uniform_real_distribution<double>(medium_uniform_left, medium_uniform_right);
			break;
		case HIGH_QUALITY:
			uniform_dist = std::uniform_real_distribution<double>(high_uniform_left, high_uniform_right);
			break;
		default:
			break;
		}

		stddev = abs(uniform_dist(generator_main));
	}

	// for annotator k, let him repeatly vote for the pair-wise comparison;
	// i.e., (o_i, o_j) , due to the wrong answer probability, especially 
	// for similar objects.
	for (int r = 0; r < num_repeat; ++r){
		// low rank, means high score, thus it will be preferred.
		bool pref;
		double a = this->params.alpha;
		double b = this->params.beta;

		if (!Is_Beta_Simu_data_type)
			pref = generate_sim_vote_for_one_task(
			n, rank_i, rank_j, stddev);
		else
			pref = generate_sim_vote_for_one_task_for_CrowdBT(
			n, rank_i, rank_j, a, b);
		if (!pref){ //  (o_i < o_j);
				if (IsDisplay){
					cout << " ----  " << r << " voting : " << " o_" << j
						<< " > o_" << i << ", ";
				}
				// make sure "i" means the better one;
				std::swap(i, j);
			}
		else{
				if (IsDisplay){
					cout << " ----  " << r << " voting : " << " o_" << i
						<< " > o_" << j << ", ";
				}
			}
		
		//**********************************
		// once the worker finishes voting to generate a preference result,
		// for example (o_j > o_i). Then based on this new observation,
		// we update the parameters according to
		// E.q. (12),(13),(14),(15),(17) and (18).
		//**********************************
		double c1 = 0, c = 0;
		NEW_PARAMS new_para = this->get_updated_parameters(i, j, c1, c);

		// update parameters 
		this-> keep_new_params(new_para, i, j);
	} /*end of repeated voting for one worker.*/

}


	// calculate the Kendall-tau distance between two ranking:
	//  1) ranking 1: is the input, usually the ground truth ranking or
	//     the result obtained from other algorithms;
	//  2) ranking 2: is the ranking obtained by this algorithm. 
double Crowd_BT::get_kendall_tau_distance(
		// ranking 1
		const int * p_rank_gt, /*ground truth ranking*/
		const bool & IsDisplay){

		// ranking 2;
		int * p_rank_crowd_bt = get_ascending_ranking();
		if (IsDisplay){
			std::cout << "The ranking:\n" << "  -- Ground Truth : ";

			for (int i = 0; i < num_o - 1; ++i){
				cout << p_rank_gt[i] << " < ";
			}

			cout << p_rank_gt[num_o - 1] << endl;

			std::cout << "  -- Crowd_BT Ranking : ";
			for (int i = 0; i < num_o - 1; ++i){
				cout << p_rank_crowd_bt[i] << " < ";
			}

			cout << p_rank_crowd_bt[num_o - 1] << endl;
		}
		double tau2 = kendallSmallN(p_rank_gt, "GT", p_rank_crowd_bt, 
			"Crowd_BT", num_o, false);
		if (IsDisplay){
			cout << "  -- Kendall-Tau distance is " << tau2 << "\n\n";
		}
		delete [] p_rank_crowd_bt;
		return tau2;
	}