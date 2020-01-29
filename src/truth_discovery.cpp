/**
* @file: truth_discovery.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "truth_discovery.hpp"
using namespace td_wieghts;

void truth_discovery::free_all(){
	int i;
	for (i = 0; i < MAXCLIENT + 1; i++){
		free(N_s[i]);
		free(G_t[i]);
	}
	free(N_s);
	free(G_t);

	for (i = 0; i < MAXOBJ + 1; i++){
		free(prob_map[i]);
		free(result_map[i]);
	}
	free(prob_map);
	free(result_map);

	for (i = 0; i < MAXTASK + 1; i++){
		free(task_list[i]);
		free(S_n[i]);
	}
	free(task_list);
	free(S_n);

	free(taskn_of_client);
	free(clientn_of_task);
	free(isExpert);
}

void truth_discovery::init(const int & first_time, const int & num_object, const int & num_worker,
	const int & max_task_num){

	MAXTASK = max_task_num;
	task_n = 0; obj_n = num_object; client_n = num_worker;
	
	int i, j;

	if (first_time){
		G_t = (int**)malloc(sizeof(int*)*(MAXCLIENT + 1));
		N_s = (client_vote**)malloc(sizeof(client_vote*)*(MAXCLIENT + 1));
	}
	for (i = 0; i < MAXCLIENT + 1; i++){
		if (first_time){
			N_s[i] = (client_vote*)malloc(sizeof(client_vote)*(MAXTASK + 1));
			G_t[i] = (int*)malloc(sizeof(int)*(MAXTASK + 1));
		}
		memset(G_t[i], 0, sizeof(int)*(MAXTASK + 1));
		memset(N_s[i], 0, sizeof(client_vote)*(MAXTASK + 1));
	}

	if (first_time){
		prob_map = (float**)malloc(sizeof(float*)*(MAXOBJ + 1));
		result_map = (vote_result**)malloc(sizeof(vote_result*)*(MAXOBJ + 1));
	}
	for (i = 0; i < MAXOBJ + 1; i++){
		if (first_time){
			prob_map[i] = (float*)malloc(sizeof(float)*(MAXOBJ + 1));
			result_map[i] = (vote_result*)malloc(sizeof(vote_result)*(MAXOBJ + 1));
		}
		memset(result_map[i], 0, sizeof(vote_result)*(MAXOBJ + 1));
		memset(prob_map[i], 0, sizeof(float)*(MAXOBJ + 1));
	}

	if (first_time){
		task_list = (int**)malloc(sizeof(int*)*(MAXTASK + 1));
		S_n = (int**)malloc(sizeof(int*)*(MAXTASK + 1));
	}
	for (i = 0; i < MAXTASK + 1; i++){
		if (first_time){
			S_n[i] = (int*)malloc(sizeof(int)*(MAXCLIENT + 1));
			task_list[i] = (int*)malloc(sizeof(int) * 2);
		}
		memset(task_list[i], 0, sizeof(int) * 2);
		memset(S_n[i], 0, sizeof(int)*(MAXCLIENT + 1));
	}

	if (first_time){
		taskn_of_client = (int*)malloc(sizeof(int)*(MAXCLIENT + 1));
		clientn_of_task = (int*)malloc(sizeof(int)*(MAXTASK + 1));
		isExpert = (int*)malloc(sizeof(int)*(MAXCLIENT + 1));
		x = (float*)malloc(sizeof(float)*(MAXTASK + 1));
		w = (float*)malloc(sizeof(float)*(MAXCLIENT + 1));
	}

	memset(taskn_of_client, 0, sizeof(int)*(MAXCLIENT + 1));
	memset(clientn_of_task, 0, sizeof(int)*(MAXTASK + 1));
	memset(isExpert, 0, sizeof(int)*(MAXCLIENT + 1));
	memset(w, 0, sizeof(float)*(MAXCLIENT));
	memset(x, 0, sizeof(float)*(MAXTASK + 1));
	
}

// Bootstrapping up to B times;
// Equation (2) in this paper;
// int update_xhat_withB(const int & B){
int truth_discovery::update_xhat_withB(const int & B, const double & _ACCURCY){
	int i, j, k, l, m, tmpid, tmpwho, x_sn;
	int X_bn_c = 0, accurcy;
	float up, down;

	int *S_bn = (int*)malloc(sizeof(int)*(MAXCLIENT + 1));
	float *xhat = (float*)malloc(sizeof(float)*(MAXTASK + 1));
	float *x_boot = (float*)malloc(sizeof(float)*(MAXTASK + 1));

	for (i = 1; i <= task_n; i++){
		x_boot[i] = 0;

		up = 0; down = 0;
		for (k = 1; k <= clientn_of_task[i]; k++){
			tmpwho = G_t[S_n[i][k]][i];
			tmpid = i;
			if (tmpwho == task_list[tmpid][0])
				x_sn = 1;
			else
				x_sn = 0;
			up += 1.0*w[S_n[i][k]] * x_sn;
			down += 1.0*w[S_n[i][k]];
		}
		x_boot[i] += 1.0*up / down;
	}

	accurcy = 0;
	for (i = 1; i <= task_n; i++){
		if (abs(x_boot[i] - x[i]) > accurcy)
			accurcy = abs(x_boot[i] - x[i]);
		x[i] = x_boot[i];
	}

	free(S_bn);
	free(xhat);
	free(x_boot);

	if (accurcy < _ACCURCY)return 1;
	return 0;
}

// Equation (1) in this paper;
/*具体改动为：
* 增加了专家标志数组 isExpert[];
* 在权值更新操作update_w_withAlpha()中加入了二次Normalization
*/
int truth_discovery::update_w_withAlpha(const double & alpha, const double & _ACCURCY){
	int i, j, k, l, m, tmpwho;
	float left, right, mid, up, down, accuracy;

	float *_w = (float*)malloc(sizeof(float)*(MAXCLIENT + 1));
	int *expert_list = (int*)malloc(sizeof(int)*(MAXCLIENT + 1));

	for (i = 1; i <= client_n; i++){
		left = 0; right = 10000;
		while (right - left > 0.001){
			mid = (left + right) / 2;
			if (boost::math::gamma_p(0.5*taskn_of_client[i], 0.5*mid) > (1 - alpha / 2))
				right = mid;
			else
				left = mid;
		}
		up = mid; down = 0;
		for (j = 1; j <= taskn_of_client[i]; j++){
			if (N_s[i][j].who == task_list[N_s[i][j].task_id][0])
				tmpwho = 1;
			else
				tmpwho = 0;
#if 0
			down += (tmpwho - x[N_s[i][j].task_id])*(tmpwho - x[N_s[i][j].task_id]);
#else
			down += abs(tmpwho - x[N_s[i][j].task_id]);
#endif
		}
		if (down == 0){
			isExpert[i] = 1;
			_w[i] = 1;
			continue;
		}
		_w[i] = up / down;
	}

	float total = 0;
	for (i = 1; i <= client_n; i++){
		if (isExpert[i] == 0)
			total += _w[i];
	}

	for (i = 1; i <= client_n; i++){
		if (isExpert[i])continue;
		_w[i] /= total;
	}

	total = 0;
	for (i = 1; i <= client_n; i++)
		total += _w[i];

	accuracy = 0;
	for (i = 1; i <= client_n; i++){
		_w[i] /= total;
		if (abs(w[i] - _w[i]) > accuracy){
			accuracy = abs(w[i] - _w[i]);
		}
		w[i] = _w[i];
	}

	free(_w);
	free(expert_list);

	if (accuracy < _ACCURCY)return 1;
	return 0;
}


void truth_discovery::iteration_main( const double & _ACCURCY,
	const int & totalTime, int & count)
{
	int i, r1, r2;
	for (i = 1; i <= client_n; i++)
		w[i] = 1.0f / client_n;

	count = 0;
	while (count < totalTime){
		count++;
		r2 = update_xhat_withB(100, _ACCURCY);
		r1 = update_w_withAlpha(0.1, _ACCURCY);
		if (r1 && r2) break;
		//std::cout << " text\n";
		//if (r1 || r2)break;
	}
}


/* Inputs:
 *   obj1: the index of object 1;
 *   obj2: the index of object 2;
 * This function is used to store the tasks, each of which consists of
 * a pair of objects for preference voting.
 * 必须保证物体的编号按照顺序来（即1~10这样连续的号）。
 * 不必输入task 的数量，每次存入会自动增加task的数量。
 * 将按照存入顺序分配task_id。
 */
void truth_discovery::store_task_pair(const int & obj1, const int & obj2){
	task_n++;
	task_list[task_n][0] = obj1;
	task_list[task_n][1] = obj2;
}

/*
 * Inputs:
 *   - worker：the index of worker;
 *   - task：the index of task.
 *   - decision：投了谁，必须保证decision是对应task的obj1或obj2中的一个。
 * 此接口用于存入用户的决策信息，同样不必定义task和worker的数量，
 * 但要保证worker的编号连续性。
 */
void truth_discovery::store_worker_decisions(int worker, int task, int decision){

	taskn_of_client[worker]++;
	N_s[worker][taskn_of_client[worker]].task_id = task;
	N_s[worker][taskn_of_client[worker]].who = decision;
	clientn_of_task[task]++;
	S_n[task][clientn_of_task[task]] = worker;
	// newly added. This is the bug, finally! about iteration times.
	G_t[worker][task] = decision;
}

void truth_discovery::print_worker_decisions(std::ofstream & of){
	of << client_n << " " << obj_n << "\n";
	for (int k = 1; k <= client_n; k++){
		of << taskn_of_client[k] << " ";
		for (int i = 1; i <= taskn_of_client[k]; i++){
			int t_id = N_s[k][i].task_id;
			std::cout << t_id << ": " << task_list[t_id][0] << ", "
				<< task_list[t_id][1] << ", " << N_s[k][i].who << ";\t ";
			of << t_id << " " << N_s[k][i].who << " ";
		}
		std::cout << std::endl;
		of << std::endl;
	}

	of << task_n << std::endl;
	for (int i = 1; i <= task_n; ++i)
		of << task_list[i][0] << " "
		<< task_list[i][1] << "\n";
}

float * truth_discovery::load_w(const bool & isDisplay){
	if (isDisplay){
		std::cout << "Weights of workers:\n";
		for (i = 1; i <= client_n; i++)
			printf("%f\t", w[i]);
		printf("\n");
	}
	return w;
}
float * truth_discovery::load_x(const bool & isDisplay){
	if (isDisplay){
		std::cout << "Scores of pairwise comparison tasks:\n";
		for (i = 1; i <= task_n; i++)
			printf("%f\t", x[i]);
		printf("\n");
	}
	return x;
}