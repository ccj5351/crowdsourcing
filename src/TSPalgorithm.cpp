/**
* @file: TSPalgorithm.cpp
* @brief: implement the SA(simulation annealing algorithm to find optimum
Hamiltonian paths.)
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 15-04-2016
*/
#include "TSPalgorithm.hpp"
#include <algorithm>

const int rand_value = 100;

TSPalgorithm::TSPalgorithm(const int & n, double ** mat, const int & idx) : permuStartingIdx(idx)
{
	tour = new Tour(n, mat);
}


TSPalgorithm::~TSPalgorithm(void)
{
}


// Set up TSP problem to run
void TSPalgorithm::Initialize(const int & iter,
	const float& temp,
	const float& cool)
{

	iterations = iter;
	temperature = temp;
	cooling_rate = cool;
	AddObserver(this);
}

// set weight matrix; 
void TSPalgorithm::setMatrix(const int & n, double ** mat){
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			tour->matrix[i][j] = mat[i][j];
}


// Run the optimization algorithm
std::vector<int> TSPalgorithm::Run(const int & flag,
	double & score, const int &startingCity,
	const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
	)
{
	tour->ResetTour();
	if (flag == 0)
		tour->CreateRandomTour();
	else if (flag == 1)
		tour->CreateNearestNeighbourTour(startingCity);
	else if (flag == 2)
		tour->RankBasedonWeights(startingCity);
	else if (flag == 3)
		tour->RankBasedonWeights();
	else if (flag == 4)
		tour->CondorcetFusionTour(startingCity, v_M);
	else
		tour->CreateNearestNeighbourTour(startingCity);

	size = tour->TourSize();
	/* initialize random seed: */
	/* since we will use the rand() function*/
	srand(time(NULL));
	// Create the optimal tour
	double best = std::numeric_limits<double>::max();
	std::vector<int> v_temp;
	for (int i = 0; i < iterations; i++)
	{
		Rotate(i);
		Reverse(i);
		SwapRandomPair(i);

		double dist = tour->HPDistance();
		// save the path;
		// the key:
		string temp;
		if (v_path.size() < TOP_K){
			for (auto i : tour->v_cities)
				temp += to_string(i);
			if (dist < LOG_ZERO_PROBABILITY){
				// the path
				v_path[temp] = tour->v_cities;
			}
		}

		if (dist < best)
		{
			best = dist;
			Notify(tour->HPDistance(), (const int&)i,
				v_temp);
		}
	}
	// cout << " ******* The best is:\n";
	Notify(best, iterations - 1, v_temp);
	score = best;
	return v_temp;
}


// Run the optimization algorithm
std::vector<int> TSPalgorithm::Run(const int & flag, const int &startingCity,
double & ResultScore, bool & ResultIsHP,
double & initialScore, bool & initialIsHP,
const std::vector<std::vector<std::vector<int>>> & v_M // each voter's voting result;
){
	tour->ResetTour();
	if (flag == 0)
		tour->CreateRandomTour();
	else if (flag == 1)
		tour->CreateNearestNeighbourTour(startingCity);
	else if (flag == 2)
		tour->RankBasedonWeights(startingCity);
	else if (flag == 3)
		tour->RankBasedonWeights();
	else if (flag == 4)
		tour->CondorcetFusionTour(startingCity, v_M);
	else
		tour->CreateNearestNeighbourTour(startingCity);
	initialScore = tour->HPDistance();
	initialIsHP = tour->isHP();
	size = tour->TourSize();
	/* initialize random seed: */
	/* since we will use the rand() function*/
	srand(time(NULL));
	// Create the optimal tour
	double best = std::numeric_limits<double>::max();
	std::vector<int> v_temp;
	for (int i = 0; i < iterations; i++)
	{
		Rotate(i);
		Reverse(i);
		SwapRandomPair(i);

		double dist = tour->HPDistance();
		// save the path;
		// the key:
		string temp;
		if (v_path.size() < TOP_K){
			for (auto i : tour->v_cities)
				temp += to_string(i);
			if (dist < LOG_ZERO_PROBABILITY){
				// the path
				v_path[temp] = tour->v_cities;
			}
		}

		if (dist < best)
		{
			best = dist;
			Notify(tour->HPDistance(), (const int&)i,
				v_temp);
		}
	}
	// cout << " ******* The best is:\n";
	Notify(best, iterations - 1, v_temp);
	ResultScore = best;
	ResultIsHP = tour -> isHP(v_temp);
	return v_temp;
}

// Get the tour node, given the index
int TSPalgorithm::GetTourCity(const int& index)
{
	int node = tour->GetCity(index);
	return node;
}


// Notify observers of changes
void TSPalgorithm::Notify(const double& dist, const int& iter,
	std::vector<int> & v)
{
	for (std::vector<Observer*>::iterator it = ob_set.begin();
		it != ob_set.end();
		it++){
		Observer* ob = *it;
		tour->Copy(v);
		//ob->UpdateDisplay(dist, iter, v);
	}
}


// Add observer
void TSPalgorithm::AddObserver(Observer* ob)
{
	ob_set.push_back(ob);
}

// Do 2-opt combinations
void TSPalgorithm::TwoOpt()
{
	double dist = tour->HPDistance();

	// Get feasible end points
	int c1_min = 1;
	int c1_max = size - 5;
	int c1;

	do
	{
		c1 = rand() % (size);
	} while (c1 < c1_min || c1 > c1_max);

	int c2 = c1 + 1;
	int c3_min = c1 + 3;
	int c3 = c3_min + rand() % (size - 1 - c3_min);
	int c4 = c3 + 1;

	// Store the city values
	int c1s = tour->GetCity(c1);
	int c2s = tour->GetCity(c2);
	int c3s = tour->GetCity(c3);
	int c4s = tour->GetCity(c4);

	tour->DoTwoOpt(c1, c2, c3, c4);

	double new_dist = tour->HPDistance();

	bool accept = Accept(new_dist, dist);

	if (!accept)
	{
		//tour.SetCities( store );	
		tour->SetCity(c1, c1s);
		tour->SetCity(c2, c2s);
		tour->SetCity(c3, c3s);
		tour->SetCity(c4, c4s);

		dist = tour->HPDistance();
	}
}

// Do 3-opt combinations
void TSPalgorithm::ThreeOpt(const int& iter)
{
	// Get tour size
	int size = tour->TourSize();

	for (int i = 0; i < size - 5; i++)
	{
		for (int j = i + 3; j < size - 3; j++)
		{
			for (int k = j + 3; k < size - 1; k++)
			{
				int c1 = i;
				int c2 = i + 1;
				int c3 = j;
				int c4 = j + 1;
				int c5 = k;
				int c6 = k + 1;

				tour->DoThreeOpt(c1, c2, c3, c4, c5, c6);
			}
		}
		std::vector<int> v_temp;
		Notify(tour->HPDistance(),
			(const int&)iter, v_temp);
	}
}


// Perform rearrangement of chosen cities:
// choose section of path and reverse the order
void TSPalgorithm::Reverse(const int& iter)
{
	// Store the tour	
	std::vector<int> store;
	tour->Copy(store);

	double dist = tour->HPDistance();

	//int c1 = rand() % (size);
	int c1 = permuStartingIdx + rand() % (size - permuStartingIdx);

	int c2;

	do
	{
		//c2 = rand() % (size);
		c2 = permuStartingIdx + rand() % (size - permuStartingIdx);
	} while (c1 == c2);

	tour->ReverseCities(c1, c2);

	double new_dist = tour->HPDistance();

	bool accept = Accept(new_dist, dist);

	if (!accept)
	{
		tour->SetCities(store);
	}
}


// Perform rearrangement of chosen cities:
// choose sections of path and rotate the order
void TSPalgorithm::Rotate(const int& iter)
{
	// Store the tour	
	std::vector<int> store;
	tour->Copy(store);

	double dist = tour->HPDistance();

	int c[3];

	c[0] = permuStartingIdx + rand() % (size - permuStartingIdx);
	c[1] = permuStartingIdx + rand() % (size - permuStartingIdx);
	c[2] = permuStartingIdx + rand() % (size - permuStartingIdx);

	while (c[0] == c[1] ||
		c[0] == c[2] ||
		c[1] == c[2])
	{
		c[0] = rand() % (size);
		c[1] = rand() % (size);
		c[2] = rand() % (size);
	}

	std::sort(c, c + 3);
	tour->RotateCities(c[0], c[1], c[2]);

	double new_dist = tour->HPDistance();

	bool accept = Accept(new_dist, dist);

	if (!accept)
	{
		tour->SetCities(store);
	}
	else
	{
		if (new_dist < dist)
		{
			std::vector<int> v_temp;
			Notify((const int&)new_dist, (const int&)iter, v_temp);
		}
	}
}


// Do pairwise swap operator on random selected city pairs
void TSPalgorithm::SwapRandomPair(const int& iter)
{
	int c1 = permuStartingIdx + rand() % (size - permuStartingIdx);
	int c2;

	do
	{
		c2 = permuStartingIdx + rand() % (size - permuStartingIdx);
	} while (c1 == c2);

	double dist = tour->HPDistance();

	int tmp1 = tour->GetCity(c1);
	int tmp2 = tour->GetCity(c2);

	tour->SetCity(c1, tmp2);
	tour->SetCity(c2, tmp1);

	double new_dist = tour->HPDistance();

	bool accept = Accept(new_dist, dist);

	if (!accept)
	{
		tour->SetCity(c1, tmp1);
		tour->SetCity(c2, tmp2);
	}
	else
	{
		if (new_dist < dist)
		{
			std::vector<int> v_temp;
			Notify((const int&)new_dist, (const int&)iter,
				v_temp);
		}
	}
}


// Do pairwise swap operator on all selected city pairs
void TSPalgorithm::RepairAllPairs()
{
	for (int i = 1; i < size - 4; i++)
	{
		for (int j = i + 3; j < size - 1; j++)
		{
			int c1 = i;
			int c2 = i + 1;
			int c3 = j;
			int c4 = j + 1;

			double dist = tour->HPDistance();

			int tmp1 = tour->GetCity(c1);
			int tmp2 = tour->GetCity(c2);
			int tmp3 = tour->GetCity(c3);
			int tmp4 = tour->GetCity(c4);

			tour->DoTwoOpt(c1, c2, c3, c4);

			double new_dist = tour->HPDistance();

			bool accept = Accept(new_dist, dist);

			if (!accept)
			{
				tour->SetCity(c1, tmp1);
				tour->SetCity(c2, tmp2);
				tour->SetCity(c3, tmp3);
				tour->SetCity(c4, tmp4);
			}
			else
			{
				int p = 0;
			}
		}
	}

}



// Returns true if we should accept the new solution,
// false otherwise.  Use standard thermodynamics equation
// P = exp( -(old-new) / kt ) where k = Boltzmann's constant
bool TSPalgorithm::Accept(const double& new_dist,
	const double& old_dist)
{
	// Calculate acceptance probability: declines with time
	double diff = (double)-1.0 * (new_dist - old_dist);
	temperature = temperature * cooling_rate;

	// Calculate probability of accepting worse solutions
	// This probability declines with time.
	double P = exp(diff / temperature);

	if (P > 1.0) return true;
	else{
		// Generate random number between 0 and 1
		double tem_p = rand() / double(RAND_MAX);
		return tem_p < P;
	}
}

void TSPalgorithm::printHP(){
	if (v_path.empty())
		cout << "No HP found.\n";
	else {
		cout << "Possible HPs include:\n";

		for (auto i : v_path){
			cout << i.first << " : ";
			if (tour->isHP(i.second))
				cout << " is an HP, ";
			else
				cout << " is not an HP, ";
			for (auto j : i.second)
				cout << j << ", ";
			cout << ", its score = "
				<< tour->HPDistance(i.second) << endl;
		}
	}
}

void TSPalgorithm::printHP(std::ofstream & of){
	if (v_path.empty()) 
		of << "\nNo HP found.\n";
	else {
		of << "\nHP candidates include:\n";
		for (auto i : v_path){
			of << i.first << " : ";
			if (tour->isHP(i.second, of))
				of << " is an HP, ";
			else
				of << " is not an HP, ";
			for (auto j : i.second)
				of << j << ", ";
			of << ", its score = "
				<< tour->HPDistance(i.second) << endl;
		}
	}
}

void Observer::UpdateDisplay(const double & dist,
	const int & iter, const std::vector<int> & path)
{
	std::cout << "Iteration: " << iter + 1 << ", "
		<< "Total distance: " << dist << "\n";
	std::cout << "Path is : ";
	for (auto i : path)
		std::cout << i << ", ";
	std::cout << std::endl;
}