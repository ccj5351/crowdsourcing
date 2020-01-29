/**
* @file: tour.hpp
* @brief: a path or a tour, used by SA(simulation annealing algorithm).
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#include "Tour.hpp"
#include <algorithm>
#include <functional>
#include <vector>
#include <utility>/*swap()*/
#include <iterator>
#include "Baseline2.hpp"/*CondorcetFuseSort()*/
using namespace std;

Tour::~Tour(void){}
Tour::Tour(const int & n, double ** mat) : vertexNum(n)
{
	matrix = mat;
}

// Do standard 2-opt operation between selected cities
void Tour::DoTwoOpt(const int& c1,
	const int& c2,
	const int& c3,
	const int& c4)
{
	// Feasible exchanges only
	if (c3 == c1 || c3 == c2 || c4 == c1 || c4 == c2) return;

	// Leave values at positions c1, c4
	// Swap c2, c3 values
	// c3 -> c2
	// c2 -> c3
	int tmp = v_cities[c2];
	v_cities[c2] = v_cities[c3];
	v_cities[c3] = tmp;
}


// Do standard 3-opt operation 
// between selected cities;
void Tour::DoThreeOpt(const int& c1,
	const int& c2,
	const int& c3,
	const int& c4,
	const int& c5,
	const int& c6){
	// Store original orderings
	int c1old = v_cities.at(c1);
	int c2old = v_cities.at(c2);
	int c3old = v_cities.at(c3);
	int c4old = v_cities.at(c4);
	int c5old = v_cities.at(c5);
	int c6old = v_cities.at(c6);

	double old_distance = HPDistance();
	// Swap 1: 40526
	// c2 -> c4	
	// c4 -> c6
	// c6 -> c2
	int tmp1 = v_cities.at(c2);
	v_cities[c2] = v_cities.at(c4);
	v_cities[c4] = v_cities.at(c6);
	v_cities[c6] = tmp1;

	double new_distance1 = HPDistance();
	// Re-instate if not an improvement
	if (new_distance1 > old_distance){
		v_cities[c1] = c1old;
		v_cities[c2] = c2old;
		v_cities[c3] = c3old;
		v_cities[c4] = c4old;
		v_cities[c5] = c5old;
		v_cities[c6] = c6old;
	}
	else{
		c1old = v_cities.at(c1);
		c2old = v_cities.at(c2);
		c3old = v_cities.at(c3);
		c4old = v_cities.at(c4);
		c5old = v_cities.at(c5);
		c6old = v_cities.at(c6);
		old_distance = new_distance1;
	}


	// Swap 2: 
	// c2 -> c4	
	// c4 -> c5
	// c5 -> c2
	tmp1 = v_cities.at(c2);
	v_cities[c2] = v_cities.at(c4);
	v_cities[c4] = v_cities.at(c5);
	v_cities[c5] = tmp1;

	new_distance1 = HPDistance();

	// Re-instate if not an improvement;
	if (new_distance1 > old_distance){
		v_cities[c1] = c1old;
		v_cities[c2] = c2old;
		v_cities[c3] = c3old;
		v_cities[c4] = c4old;
		v_cities[c5] = c5old;
		v_cities[c6] = c6old;
	}
}

// Get Hamiltonian Path Distance;
// HP distance := sigma(log10(1/w_ij));
// where, w_ij is the weight of directed edge (i, j),
// i.e., the edge (i --> j);
double Tour::TourDistance()const{
	double dist = 0.0;
	int size = (int)v_cities.size();
	int c1, c2;
	for (int i = 0; i < size - 1; i++)
	{
		c1 = v_cities.at(i);
		c2 = v_cities.at(i + 1);
		dist += Log10Distance(c1, c2);
	}
	// And back to the beginning city;
	// that is, the Hamiltonian Cycle;
	c1 = v_cities.at(size - 1);
	c2 = v_cities.at(0);
	dist += Log10Distance(c1, c2);
	return dist;
}

// Get total distance of tour, 
// that is the Hamiltonian Cycle;
// distance := sigma(log10(1/w_ij)).
/*double Tour::HPDistance()const {
double dist = 0.0;
int size = (int)v_cities.size();
int c1, c2;
for (int i = 0; i < size - 1; i++)
{
c1 = v_cities.at(i);
c2 = v_cities.at(i + 1);
if (0.0 != matrix[c1][c2])
continue;
else
dist += 1.0;
}
return dist;
}*/

bool Tour::isHP(const std::vector<int> & v)const {
	int size = (int)v.size();
	for (int i = 0; i < size - 1; i++)
	{
		if (matrix[v[i]][v[i + 1]] == 0){
			/*std::cout << " (" << v[i] << ", " << v[i + 1] <<
				") = 0, ";*/
			return false;
		}
	}
	return true;
}

bool Tour::isHP()const {
	int size = (int)this->v_cities.size();
	for (int i = 0; i < size - 1; i++)
	{
		if (this->matrix[v_cities[i]][v_cities[i + 1]] == 0){
			/*
			std::cout << " (" << v_cities[i] << ", " << v_cities[i + 1] <<
				") = 0, ";*/
			return false;
		}
	}
	return true;
}

bool Tour::isHP(const std::vector<int> & v, std::ofstream & of)const {
	int size = (int)v.size();
	for (int i = 0; i < size - 1; i++)
	{
		if (matrix[v[i]][v[i + 1]] == 0){
			/*
			of << " edge (" << v[i] << ", " << v[i + 1] <<
				") = 0, ";*/
			return false;
		}
	}
	return true;
}


double Tour::HPDistance()const {
	double dist = 0.0;
	int size = (int)v_cities.size();
	int c1, c2;
	for (int i = 0; i < size - 1; i++)
	{
		c1 = v_cities.at(i);
		c2 = v_cities.at(i + 1);
		dist += Log10Distance(c1, c2);
	}
	return dist;
}

double Tour::HPDistance(const std::vector<int> & v) const{
	double dist = 0.0;
	int size = (int)v.size();
	int c1, c2;
	for (int i = 0; i < size - 1; i++)
	{
		c1 = v.at(i);
		c2 = v.at(i + 1);
		dist += Log10Distance(c1, c2);
	}
	return dist;
}

// Generate arbitrary tour;
void Tour::CreateRandomTour(){
	ResetTour();

	for (int i = 0; i < vertexNum; i++){
		v_cities.push_back(i);
	}
	// keep the first vertex unchanged;
	// thus, using "begin() + 1" instead of "begin()".
	// random_shuffle( v_cities.begin() + 1, v_cities.end() );
	random_shuffle(v_cities.begin(), v_cities.end());
}

// Generate arbitrary tour;
void Tour::CreateTour(){
	ResetTour();
	for (int i = 0; i < vertexNum; i++){
		v_cities.push_back(i);
	}
}

// Get distance between selected cities;
double Tour::Log10Distance(const int& c1,
	const int& c2) const
{
	double w = matrix[c1][c2] + DELTA_PROBABILITY;
	return (-1.0)*log10(w);
}

// Set pointer to cost / distance matrix object;
void Tour::SetMatrix(double ** mat){
	matrix = mat;
}

// Reset existing tour data;
void Tour::ResetTour(){
	std::vector<int>().swap(v_cities);
}

// Return the number of cities in the tour;
int Tour::TourSize() const{
	return vertexNum;
}

// Get the tour city, given the index passed to it;
int Tour::GetCity(const int& index){
	return v_cities[index];
}

// Get tour from the set of nearest neighbors;
void Tour::CreateNearestNeighbourTour(const int & startingCity){
	ResetTour();
	std::set<int> city_set;
	std::set<int>::iterator it;

	for (int i = 0; i < vertexNum; i++)
	{
		city_set.insert(i);
	}

	int city = startingCity;

	for (int i = 1; i <= vertexNum; i++)
	{
		// Add city to tour
		v_cities.push_back(city);

		// Remove city from unique set
		it = city_set.find(city);
		city_set.erase(it);

		city = GetNearestNeighbour(city, city_set);
	}
}

// similar to weighted CondorcetFuseTour,
// but this method may not be correct,
// due to overfitting problem.
void Tour::RankBasedonWeights(){
	//calculate the weight of each node
	vector<std::pair<int, double>> v_weight(vertexNum);
	for (int i = 0; i < vertexNum; ++i){
		v_weight[i] = make_pair(i, .0);
	}

	// in ----> +;
	// out ---> -;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_weight[i].second += matrix[j][i];
		}
	}

	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_weight[i].second -= matrix[i][j];
		}
	}

	std::sort(v_weight.begin(), v_weight.end(),
		[](const std::pair<int, double> & arg1,
		const std::pair<int, double> & arg2){return arg1.second < arg2.second; });

	if (v_cities.empty()){
		for (int i = 0; i < vertexNum; ++i)
			v_cities.push_back(i);
	}
	for (int i = 0; i < vertexNum; ++i){
		v_cities[i] = v_weight[i].first;
	}

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	cout << "Baseline ranking result is :\n";
	for (int i = 0; i < vertexNum; ++i){
		cout << p_rank[i] << ",  ";
	}
	cout << std::endl;
#endif

	vector<std::pair<int, double>>().swap(v_weight);
}


void Tour::CondorcetFusionTour(const int & startingCity, 
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	){
	// each of the w search systems, or w voters; 
	const int w = v_M.size();
	int * arr = new int[this->vertexNum];
	for (int i = 0; i < vertexNum; ++i)
		arr[i] = i;
	CondorcetFuseSort(arr, 0, vertexNum - 1, vertexNum, w, v_M );
	if (this -> v_cities.empty()){
		for (int i = 0; i < vertexNum; ++i)
				this-> v_cities.push_back(i);
		}

	for (int i = 0; i < vertexNum; ++i){
		// reverse the tour, since we have different ranking strategy,
		// i.e., for the simulated data, i --> j, means prefer i than j;
		// which is used for RC, GM.

		// but I used to the notation , i--->j, ,menas prefer j than i.
		// which is used for SA.
		v_cities[i] = arr[vertexNum - 1 - i];
	}
	
	auto it = find(v_cities.begin(), v_cities.end(), startingCity);
	auto idx = std::distance(v_cities.begin(), it);
	// swap two vertices;
	v_cities[idx] = v_cities[0];
	v_cities[0] = startingCity;
	delete [] arr;
}

void Tour::RankBasedonWeights(const int & startingCity){
	//calculate the weight of each node
	vector<std::pair<int, double>> v_weight(vertexNum);
	for (int i = 0; i < vertexNum; ++i){
		v_weight[i] = make_pair(i, .0);
	}

	// in ----> +;
	// out ---> -;
	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_weight[i].second += matrix[j][i];
		}
	}

	for (int i = 0; i < vertexNum; ++i){
		for (int j = 0; j < vertexNum; ++j){
			v_weight[i].second -= matrix[i][j];
		}
	}

	std::sort(v_weight.begin(), v_weight.end(),
		[](const std::pair<int, double> & arg1,
		const std::pair<int, double> & arg2){return arg1.second < arg2.second; });
	
	if (v_cities.empty()){
		for (int i = 0; i < vertexNum; ++i)
			v_cities.push_back(i);
	}

	for (int i = 0; i < vertexNum; ++i){
		v_cities[i] = v_weight[i].first;
	}

#define __RELEASE_MODE_CORWDSOURCING_PROJECT_
#ifndef __RELEASE_MODE_CORWDSOURCING_PROJECT_
	cout << "Baseline ranking result is :\n";
	for (int i = 0; i < vertexNum; ++i){
		cout << p_rank[i] << ",  ";
	}
	cout << std::endl;
#endif
	auto it = find(v_cities.begin(), v_cities.end(),startingCity);
	auto idx = std::distance(v_cities.begin(), it);
	// swap two vertices;
	v_cities[idx] = v_cities[0];
	v_cities[0] = startingCity;
	vector<std::pair<int, double>>().swap(v_weight);
}


// Get the nearest node to this one;
int Tour::GetNearestNeighbour(const int& c, std::set<int>& cset)
{
	int city = 0;
	// Get minimum distance node
	double min_dist = 99999999;
	std::set<int>::iterator cit;

	for (cit = cset.begin(); cit != cset.end(); cit++){
		int n = *cit;

		if (n == c) continue;
		double dist = Log10Distance(c, n);
		if (dist < min_dist){
			city = n;
			min_dist = dist;
		}
	}
	return city;
}

// Set the cities;
void Tour::SetCities(const std::vector<int>& v){
	v_cities = v;
}

// Set city;
void Tour::SetCity(const int& index,
	const int& value){
	v_cities[index] = value;
}


// Copy / clone the tour cities to another vector;
void Tour::Copy(std::vector<int>& v){
	v.reserve(vertexNum);
	v = v_cities;
}

// Reverse the path between the selected end points;
void Tour::ReverseCities(const int& c1, const int& c2){
	int c1_n = c1 < c2 ? c1 : c2;
	int c2_n = c1 < c2 ? c2 : c1;

	if (c1_n == 0) c1_n = 1;
	std::vector<int>::iterator it1, it2;

	it1 = v_cities.begin() + c1_n;
	it2 = v_cities.begin() + c2_n;
	reverse(it1, it2);
}

// Rotate path cities given the beginning, 
// middle and end points;
void Tour::RotateCities(const int& c1,
	const int& c2,
	const int& c3)
{
	int c1_new = c1 == 0 ? 1 : c1;

	std::vector<int>::iterator it1 = v_cities.begin() + c1_new;
	std::vector<int>::iterator it2 = v_cities.begin() + c2;
	std::vector<int>::iterator it3 = v_cities.begin() + c3;
	rotate(it1, it2, it3);
}