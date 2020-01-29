/**
* @file: TSPalgorithm.hpp
* @brief: implement the SA(simulation annealing algorithm to find optimum
          Hamiltonian paths.)
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 15-04-2016
*/
#pragma once
#include "Tour.hpp"
#include <string>
#include <iostream>
#include <fstream>      // std::ofstream
#include <map>
#include <ctime> /*rand()*/
#include <utility> // std::pair, std::make_pair
using namespace std;
#define TOP_K 10
class Observer
{
public:
	void UpdateDisplay(const double& dist,
		const int & iter, const std::vector<int> & path);
};

class TSPalgorithm : public Observer
{
public:
	TSPalgorithm(const int & n, double ** mat, const int & idx = 1);
	~TSPalgorithm(void);

	void Initialize(const int& iterations,
		const float& temp,
		const float& cool);

	// Run the optimization algorithm
	std::vector<int> Run(const int & flag, double & score, 
		const int &startingCity,
		const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
		);
	

	// Run the optimization algorithm
	std::vector<int> Run(const int & flag, const int &startingCity,
		double & score, bool & ResultIsHP,
		double & iniScore, bool & iniIsHP,
		const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
		);

	int  GetTourCity(const int& index);
	void setMatrix(const int & n, double ** mat);
	void AddObserver(Observer* ob);
	void Notify(const double& dist, const int& iter,
		std::vector<int> & v);

	bool Accept(const double& new_dist,
		const double& old_dist);
	void printHP();
	void printHP(std::ofstream & of);

private:
	const int permuStartingIdx;
	void Rotate(const int& iter);
	void Reverse(const int& iter);
	void TwoOpt();
	void SwapRandomPair(const int& iter);
	void RepairAllPairs();
	void ThreeOpt(const int& iter);
public:
	Tour * tour;
private:
	int iterations;
	float temperature;
	float cooling_rate;
	std::vector<Observer*> ob_set;
	int size;
	map<string, vector<int>> v_path;
public:
	pair<vector<int>, double> opt_path;
};

