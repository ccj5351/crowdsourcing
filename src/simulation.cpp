/**
* @file: simulation.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#include "simulation.hpp"
using namespace std;

//read the edge information matrix
bool scanMatrix(double **  matrix, const int & n,
	const string & filename){
	//scan the file
	ifstream myfile(filename.c_str());
	if (!myfile){
		cout << "failed open file " << filename
			<< endl;
		return false;
	}
	char c;
	int row = 0;
	int col = 0;
	int count = -1;
	string temp = "";
	while (myfile.get(c)){
		temp += c;
		if (c == ','){
			count++;
			//cout << temp << endl;
			for (int i = 0; i < temp.length(); ++i){
				if (temp[i] == '0'){
					matrix[row][col] = 0;
					col++;
				}
				else if (temp[i] == '1'){
					matrix[row][col] = 1;
					col++;
				}
			}
			if (count % n == (n - 1)){
				row++;
				col = 0;
			}
			temp = "";
		}
	}
	//cout << count << endl;
	myfile.close();
	return true;
}
bool scanMatrix2(double **  matrix, const int & n,
	const string & filename){
	//scan the file
	ifstream myfile(filename.c_str());
	if (!myfile){
		cout << "failed open file " << filename
			<< endl;
		return false;
	}
	char buff[100];
	char c;
	int row = 0;
	int col = 0;
	int count = -1;
	string temp = "";
	while (myfile.get(c)){
		temp += c;
		if (c == ',' || c == ';'){
			count++;
			//cout << temp << endl;
			int loc = 0;
			for (int i = 0; i < temp.length(); ++i){
				if ((temp[i] >= '0' && temp[i] <= '9') || temp[i] == '.'){
					//matrix[row][col] = 0;
					buff[loc++] = temp[i];
				}
			}

			buff[loc] = '\0';
			double db = atof(buff);
			matrix[row][col] = db;
			col++;

			if (count % n == (n - 1)){
				row++;
				col = 0;
			}
			temp = "";
		}
	}
	//cout << count << endl;
	myfile.close();
	return true;
}

//generate uniform number
void generate_uniform_num(int n, double* num){
	for (int i = 0; i < n; ++i){
		// generate a uniform random value [0, 1]
		num[i] = rand() / double(RAND_MAX);
	}
}

//uniform distribution simulation
void simulate_order(
	double ** task,
	double ** order_matrix,
	double * uniform_num,
	int vertexNum,
	int edge_num){

	int count = 0;
	int n1 = 0, n2 = edge_num, n3 = edge_num * 2;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = i + 1; j < vertexNum; ++j){
			// find an edge (u,v) need to decide directions
			// 1) only (u,v);
			// 2) only (v,u);
			// 3) both (u,v) and (v, u);
			if (task[i][j] > 0){
				count++;
				//first decide one of the three choices
				// the order is vi > vj
				if (uniform_num[n1] > 0 && uniform_num[n1] <= 0.45){
					order_matrix[i][j] = 1;
				}
				// the order is vj > vi
				else if (uniform_num[n1] > 0.45 &&  uniform_num[n1] <= 0.9){
					order_matrix[j][i] = 1;
				}

				else{
					//decide who the weight
					//a little controversial
					if (uniform_num[n3] > 0 && uniform_num[n3] <= 0.4){
						order_matrix[i][j] = 0.9;
						order_matrix[j][i] = 0.1;
					}
					//middle controversial
					else if (uniform_num[n3] > 0.4 && uniform_num[n3] <= 0.8){
						order_matrix[i][j] = 0.7;
						order_matrix[j][i] = 0.3;
					}
					//very controversial
					else{
						order_matrix[i][j] = 0.5;
						order_matrix[j][i] = 0.5;
					}
					n3++;

					//decide which vertex has more weight
					if (uniform_num[n2] < 0.5){
						double temp = order_matrix[i][j];
						order_matrix[i][j] = order_matrix[j][i];
						order_matrix[j][i] = temp;
						n2++;
					}

				}
				n1++;
			}
		}
	}
	//cout << count << endl;
}

// set edge weights as zero's.
void clear_matrix(double ** matrix, int n){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			matrix[i][j] = 0;
		}
	}
}

// to generate the simulation graph G_T;
// just run this function to get the final result.
// all the related functions defined in this file
// will be called if used.
//normal distribution simulation
void getSimuData(const int & edgeNum, // e.g. edge_num = 20;
	const int & vertexNum, // n = 10;
	const std::string & inputFileName,
	//e.g., "C:\\Users\\Mobrick\\Desktop\\simulation\\matrix.txt"
	// matrix rows are separated by ";" for MATLAB readable-data;
	// matrix rows are separated by "," for C++ readable-data;
	const std::string & outputFileM, // for store MATLAB data
	// all the data in the matrix is output 
	// as separated rows.
	const std::string & outputFileC, // for store C++ data;
	// all the data in the matrix is output 
	// as a whole line.
	const std::string & outputFileC1 // for store C++ data in one line;
	){
	double temp = 0;
	// generate 
	double * uniform_num = new double[edgeNum * 3];
	generate_uniform_num(edgeNum * 3, uniform_num);

	double **pp_g_t = new double*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		pp_g_t[i] = new double[vertexNum];
	}
	// read the data in G_T graph.
	scanMatrix(pp_g_t, vertexNum, inputFileName);

	double **order_matrix = new double*[vertexNum];
	for (int i = 0; i < vertexNum; ++i){
		order_matrix[i] = new double[vertexNum];
	}
	// set the data to zero's.
	clear_matrix(order_matrix, vertexNum);

	simulate_order(pp_g_t, order_matrix, uniform_num,
		vertexNum, edgeNum);

	ofstream outfileM(outputFileM.c_str()),
		outfileC(outputFileC.c_str()),
		outfileC1(outputFileC1.c_str());
	if (!outfileM || !outfileC || !outfileC1){
		cout << "Can not open the file to store the result.\n";
	}

	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum - 1; ++j){
			outfileM << order_matrix[i][j] << ", ";
			outfileC << order_matrix[i][j] << ", ";
			outfileC1 << order_matrix[i][j] << ", ";
		}

		outfileM << order_matrix[i][vertexNum - 1] << ";\n";
		outfileC << order_matrix[i][vertexNum - 1] << ",\n";
		outfileC1 << order_matrix[i][vertexNum - 1] << ", ";
	}
	if (outfileM.is_open())
		outfileM.close();
	if (outfileC.is_open())
		outfileC.close();
	if (outfileC1.is_open())
		outfileC1.close();
}