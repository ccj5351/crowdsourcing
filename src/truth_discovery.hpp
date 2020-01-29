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
�����������ⲿ�������ݵĽӿڣ�
һ���������ݽӿ��÷�
void store_task_pair(int obj1, int obj2)
obj1����һ������ı��;
obj2���ڶ�������ı��;

�˽ӿ����ڴ���task�ԣ���Ҫ�Ƚϵ����������壩�����б��뱣֤����ı�Ű���˳��������1~10���������ĺţ���
��������task��������ÿ�δ�����Զ�����task��������
�����մ���˳�����task_id��


void store_worker_decisions(int worker, int task, int decision)

worker�����˵�ID���
task��task��ID��ţ����յ�һ�����������
decision��Ͷ��˭�����뱣֤decision�Ƕ�Ӧtask��obj1��obj2�е�һ����

�˽ӿ����ڴ����û��ľ�����Ϣ��ͬ�����ض���task��worker����������Ҫ��֤worker�ı�������ԡ�

������ȡ���ݽӿ��÷�

float *load_x();
float *load_w();
����˼�壬���ص���w��x�����׵�ַָ�롣

����ʹ������

���ⲿ����ʱ���ȵ���init()������ʼ���ڴ档
�ڵ��� store_task_pair() �� store_worker_decisions() �ӿں����������ݡ�
�۵��� calc_initial_prob() ���������ʼ���ʾ���
�ܵ��� itertaion_main() �������е������㡣
��ʹ�� load_x() �� load_w() �õ���������
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
		 * ���뱣֤����ı�Ű���˳��������1~10���������ĺţ���
		 * ��������task ��������ÿ�δ�����Զ�����task��������
		 * �����մ���˳�����task_id��
		 */
		void store_task_pair(const int & obj1, const int & obj2);

		/*
		 * Inputs:
		 *   - worker��the index of worker;
		 *   - task��the index of task.
		 *   - decision��Ͷ��˭�����뱣֤decision�Ƕ�Ӧtask��obj1��obj2�е�һ����
		 * �˽ӿ����ڴ����û��ľ�����Ϣ��ͬ�����ض���task��worker��������
		 * ��Ҫ��֤worker�ı�������ԡ�
		 */
		void store_worker_decisions(int worker, int task, int decision);
		void print_worker_decisions(std::ofstream & of);

		float *load_w(const bool & isDisplay);
		float *load_x(const bool & isDisplay);

	}; /*end of the struct*/

} /*end of namespace definition*/
#endif