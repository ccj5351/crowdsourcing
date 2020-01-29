/**
* @file: tour.hpp
* @brief: a path or a tour, used by SA(simulation annealing algorithm).
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#pragma once
#include <vector>
#include <set>
#include <limits>
#include <iostream>
#include <fstream>

class Tour
{
#define DELTA_PROBABILITY 1.0E-300 // a very small positive floating value;
#define LOG_ZERO_PROBABILITY (-1)*log10(DELTA_PROBABILITY)
public:
	Tour(const int & n, double ** mat);
	//Tour(const int & n, double **  mat, std::vector<int> tour);
	~Tour(void);

	void DoTwoOpt(const int& c1,
		const int& c2,
		const int& c3,
		const int& c4);

	void DoThreeOpt(const int& c1,
		const int& c2,
		const int& c3,
		const int& c4,
		const int& c5,
		const int& c6);
	bool Tour::isHP() const;
	bool Tour::isHP(const std::vector<int> & v)const;
	bool Tour::isHP(const std::vector<int> & v, std::ofstream & of)const;
	double HPDistance() const;
	double HPDistance(const std::vector<int> & v) const;
	double TourDistance()const;
	int TourSize() const;
	int GetCity(const int& index);
	void SetCity(const int& index, const int& value);

	void SetMatrix(double **  mat);
	void CreateRandomTour();
	void CreateTour();
	void CreateNearestNeighbourTour(const int & startingVertex);
	void RankBasedonWeights();
	void RankBasedonWeights(const int & startingCity);
	
	void CondorcetFusionTour(const int & startingCity,
		const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
		);
	void ResetTour();
	void SetCities(const std::vector<int>& v);
	void Copy(std::vector<int>& v);
	void ReverseCities(const int& c1, const int& c2);
	void RotateCities(const int& c1,
		const int& c2,
		const int& c3);
public:
	std::vector< int > v_cities;
	double ** matrix;
private:
	double Log10Distance(const int& c1,
		const int& c2) const;
	int GetNearestNeighbour(const int& c, std::set<int>& cset);
	const int vertexNum = 5;
};

