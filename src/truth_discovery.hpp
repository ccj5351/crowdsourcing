/**
* @file: truth_discovery.hpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#ifndef __HEADER__CROWD_SOURCING_TRUTH_DISCOVERY_HPP_
#define __HEADER__CROWD_SOURCING_TRUTH_DISCOVERY_HPP_

#include <boost/math/special_functions/gamma.hpp>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>

/*
* Truth discovery, for updating weights iteratively.
* According to the paper: "Towards COnfidence in the Truth:
* A bootstrapping based Truth Discovery Approach".
* Original Code implementation is finished by Haipei Sun(esbiehs@gmail.com).
* A little update by Changjiang Cai(ccai1@stevens.edu).
*/


/* Interfaces_readme.txt:
加入了两个外部载入数据的接口：
一、存入数据接口用法
void store_task_pair(int obj1, int obj2)
obj1：第一个物体的编号;
obj2：第二个物体的编号;

此接口用于存入task对（即要比较的那两个物体），其中必须保证物体的编号按照顺序来（即1~10这样连续的号）。
不必输入task的数量，每次存入会自动增加task的数量。
将按照存入顺序分配task_id。


void store_worker_decisions(int worker, int task, int decision)

worker：工人的ID编号
task：task的ID编号，按照第一个分配的来。
decision：投了谁，必须保证decision是对应task的obj1或obj2中的一个。

此接口用于存入用户的决策信息，同样不必定义task和worker的数量，但要保证worker的编号连续性。

二、读取数据接口用法

float *load_x();
float *load_w();
顾名思义，返回的是w和x数组首地址指针。

三、使用流程

①外部调用时首先调用init()函数初始化内存。
②调用 store_task_pair() 和 store_worker_decisions() 接口函数载入数据。
③调用 calc_initial_prob() 函数计算初始概率矩阵。
④调用 itertaion_main() 函数进行迭代计算。
⑤使用 load_x() 和 load_w() 得到计算结果。
*/
namespace td_wieghts {


//#define MAXTASK	160000
#define MAXOBJ	500
#define MAXCLIENT	5000
// #define ACCURCY 1E-10
// #define TH 100000


#define DIRECT_USE_WEIGHTS  1
#define	EXPONENTIAL_USE_WEIGHTS 2
#define	LAPLACE_SMOTHING 3

	
	struct client_vote{
		int task_id, who;
		client_vote(){};
		client_vote(const int & tid, const int & whoo){
			task_id = tid;
			who = whoo;
		}
	};

	struct vote_result{
		int x_vote, y_vote;
	};

	struct truth_discovery {
		int  MAXTASK;
		FILE *fin, *fout;
		int task_n, obj_n, client_n, i;
		client_vote **N_s;
		vote_result **result_map;
		float **prob_map, **conf_matrix;
		int **G_t, **S_n, **task_list;
		float *x, *w;
		int *taskn_of_client, *clientn_of_task, *isExpert, count;

		void free_all();

		void init(const int & first_time, const int & num_object, const int & num_worker, const int & max_task_num);


		int update_xhat_withB(const int & B, const double & _ACCURCY);

		// Equation (1) in this paper;
		int update_w_withAlpha(const double & alpha, const double & _ACCURCY);


		void iteration_main(const double & _ACCURCY,
			const int & totalTime, int & count);


		/* Inputs:
		 *   obj1: the index of object 1;
		 *   obj2: the index of object 2;
		 * This function is used to store the tasks, each of which consists of
		 * a pair of objects for preference voting.
		 * 必须保证物体的编号按照顺序来（即1~10这样连续的号）。
		 * 不必输入task 的数量，每次存入会自动增加task的数量。
		 * 将按照存入顺序分配task_id。
		 */
		void store_task_pair(const int & obj1, const int & obj2);

		/*
		 * Inputs:
		 *   - worker：the index of worker;
		 *   - task：the index of task.
		 *   - decision：投了谁，必须保证decision是对应task的obj1或obj2中的一个。
		 * 此接口用于存入用户的决策信息，同样不必定义task和worker的数量，
		 * 但要保证worker的编号连续性。
		 */
		void store_worker_decisions(int worker, int task, int decision);
		void print_worker_decisions(std::ofstream & of);

		float *load_w(const bool & isDisplay);
		float *load_x(const bool & isDisplay);

	}; /*end of the struct*/

} /*end of namespace definition*/
#endif