/**
 * @file: printParams.cpp
 * @brief:
 * @author: Changjiang Cai, caicj5351@gmail.com
 * @version: 0.0.1
 * @creation date: 17-12-2015
 * @last modified: Thu 10 Mar 2016 09:32:02 AM EST
 */

#include "printParams.hpp"

using namespace std;

void print_task_assignment_parameters(
				const bool & IsSave2File,
				ofstream & myfile,
				const int & case_idx,
				const int & n,
				const int & l, 
				const vector<int> & v_m, 
				const vector<int> & v_B){
        //ofstream myfile;
	    //int t = n*(n/2 + 1);
		int t = n*n/2;
		int size_m = v_m.size();
		int size_B = v_B.size();
		double h= .0, pt = 0.0, pw = .0;

		if (!IsSave2File){
			    cout << "Case " << case_idx << " : n = " << n << ", l = " <<l << ".\n";
		        printf(" n( # of   l( # of workers)  m(# of    h(# of tasks B         pt(payment  pw(payment\n");
		        printf(" objects)  who do some task  workers)  per workers  (budget)  per task)   per worker)\n");
		        printf(" ========  ================  ========  ===========  ========  ==========  ===========\n");
		}
		else{
		        myfile << "Case " << case_idx << " : n = " 
						<< n << ", l = " <<l << ".\n"; 
				myfile << " n( # of   l( # of workers)  m(# of    h(# of tasks B         pt(payment  pw(payment\n";
				myfile << " objects)  who do some task  workers)  per workers  (budget)  per task)   per worker)\n";
				myfile << " ========  ================  ========  ===========  ========  ==========  ===========\n";
		}
		for (int i = 0 ; i < size_m ; ++ i){
		        h = double(1.0*t*l/v_m[i]);
				for(int j = 0; j < size_B ; ++ j){
				pt = (double)v_B[j]/(t*l);
			    pw = (double)v_B[j]/v_m[i];
				if (!IsSave2File){
					    printf(" %-8d ", n);
				        printf(" %-16d ", l);
				        printf(" %-8d ", v_m[i]);
				        printf(" %-11.4f ", h);
				        printf(" %-8d ", v_B[j]);
				        printf(" %-11.4f ", pt);
				        printf(" %-10.4f\n", pw);
				}
				else{
					    myfile << " %-8d "   << n;
				        myfile << " %-16dd " << l;
				        myfile << " %-8d " << v_m[i];
				        myfile << " %-11.4f " << h;
				        myfile << " %-8d " << v_B[j];
				        myfile << " %-11.4f " << pt;
				        myfile << " %-10.4f\n" << pw;
				}
				} 
				if (!IsSave2File)
						printf(" --------  ----------------  --------  -----------  --------  ----------  -----------\n");
		        else
						myfile << " --------  ----------------  --------  -----------  --------  ----------  -----------\n"; 
		}
		if (!IsSave2File)
				printf("\n\n");
		else
				myfile << "\n\n";
}

