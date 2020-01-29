#include "MyUtility.h"
using namespace std;
ullong factorial(unsigned int n){
	if (n == 0)
		return 1;
	return(n * factorial(n - 1));
}

ullong factorial(int n){
	if (n == 0)
		return 1;
	return(n * factorial(n - 1));
}

void printTiePath(const int & v, // vertex number;
	const std::vector<int> & v_tiePath, // tie-paths;
	const std::string  & fn) // output file name;
{
	int row_num = v_tiePath.size() / v;
	int ** p_Gp = new int *[row_num*v];
	for (int i = 0; i < row_num*v; ++i){
		p_Gp[i] = new int[v];
	}
	// initialize
	for (int i = 0; i < row_num*v; ++i){
		for (int j = 0; j < v; ++j){
			p_Gp[i][j] = 0;
		}
	}

	std::ofstream  of;
	of.open(fn.c_str());
	if (!of.is_open()){
		std::cout << "Error! Not open " << fn << endl;
	}
	else{
		for (int r = 0; r < row_num; ++r){
			for (int c = 0; c < v - 1; ++c){
				p_Gp[r*v + v_tiePath[r*v + c]][v_tiePath[r*v + c + 1]] = 1;
			}
		}

		for (int i = 0; i < row_num*v; ++i){
			for (int j = 0; j < v; ++j){
				std::cout << p_Gp[i][j] << ",";
				of << p_Gp[i][j] << ",";
			}
			std::cout << endl; of << endl;
		}
	}
	// release memory
	for (int i = 0; i < row_num*v; ++i){
		delete[] p_Gp[i];
	}
	delete[] p_Gp;
}

string toString(const vector<int> &p) {
	string s;

	for (vector<int>::const_iterator it = p.begin(); it != p.end(); it++) {
		s += std::to_string(*it) + " <- ";
	}

	return s;
}

void printFactorial(std::ofstream  & of, const int & n){
	ullong delta = factorial(n) / 10;
	of << "factorial(" << n << ") = " << factorial(18) << "\n"
		<< "delta = factorial(" << n << ")/10 = " << delta << "\n";
	for (int i = 0; i < 10; ++i){
		of << "[  " << i + 1 << " " << i*delta << " " << (i + 1)*delta << "  ]\n";
	}

	for (int i = 0; i < 10; ++i){
		of << i + 1 << " " << i*delta << " " << (i + 1)*delta << "\n";
	}
}


std::vector<std::vector<double> >
add(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs) {
	vector<vector<double> > result;

	result = lhs;
	for (int i = 0; i < (int)lhs.size(); i++) {
		for (int j = 0; j < (int)lhs[0].size(); j++)
			result[i][j] = lhs[i][j] + rhs[i][j];
	}

	return result;
}

std::vector<std::vector<int>> add(const std::vector<std::vector<int> >& lhs, const std::vector<std::vector<int> >& rhs){
	vector<vector<int> > result;

	result = lhs;
	for (int i = 0; i < (int)lhs.size(); i++) {
		for (int j = 0; j < (int)lhs[0].size(); j++)
			result[i][j] = lhs[i][j] + rhs[i][j];
	}
	return result;
}


std::vector<std::vector<double> >
multiply(const std::vector<std::vector<double> >& lhs, const std::vector<std::vector<double> >& rhs) {
	int dim1, dim2, dim3;
	vector<vector<double> > result;
	vector<double> row;

	//get dimensions of lhs and rhs
	dim1 = (int)lhs.size();
	dim2 = (int)lhs[0].size();
	dim3 = (int)rhs[0].size();

	//initialize result matrix
	for (int i = 0; i < dim1; i++) {
		row = vector<double>(dim3);
		result.push_back(row);
	}

	//assign value
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim3; j++) {
			result[i][j] = 0;
			for (int k = 0; k < dim2; k++)
				result[i][j] += lhs[i][k] * rhs[k][j];
		}
	}

	return result;
}

// direct weights
std::vector<std::vector<double> >&
unify_BT_D(std::vector<std::vector<double> >& M) {
	size_t dim1 = M.size();
	double sum;

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < i; j++) {
			sum = M[i][j] + M[j][i];
			if (sum != 0) {
				M[i][j] = M[i][j] / sum;
				M[j][i] = M[j][i] / sum;
			}
		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < dim1; ++i)
		M[i][i] = 0.0;

	return M;
}

// exponential weights
std::vector<std::vector<double> >&
unify_BT_E(std::vector<std::vector<double> >& M) {
	size_t dim1;
	double sum;
	dim1 = M.size();

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < i; j++) {
			sum = exp(M[i][j]) + exp(M[j][i]);
			if (sum != 0) {
				M[i][j] = exp(M[i][j]) / sum;
				M[j][i] = exp(M[j][i]) / sum;
			}
		}
	}
	// set M(i,i) = 0.0;
	for (int i = 0; i < dim1; ++i)
		M[i][i] = 0.0;

	return M;
}

void
output_matrix(std::vector<std::vector<int> > & matrix) {
	int n;

	n = (int)matrix.size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}


void
output_matrix(std::vector<std::vector<double> > &matrix) {
	int n;

	n = (int)matrix.size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}


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