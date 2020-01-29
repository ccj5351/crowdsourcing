/**
* @file: repeatchoice.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "repeatchoice.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <utility>
#include <random>
#include <algorithm>
#include <iostream>
#include <chrono>

using namespace std;

RepeatChoice::RepeatChoice(const int &_n, 
	const std::vector<std::string> &_files) {
    ifstream ifs;
    istringstream iss;
    vector<vector<int> > M;
    string line, file;
    set<pair<int, int> > partial_ranking;
    int count;
    
    this->n = _n;
    this->files = _files;
    
    for (vector<string>::iterator it = this->files.begin(); 
		it != this->files.end(); it++) {
        file = *it;
        M.clear();
        ifs.open(file);
        for (int i = 0; i < n; i++) {
            getline(ifs, line);
            std::replace(line.begin(), line.end(), ',', ' ');
            M.push_back(vector<int> (n));
            iss.str(line);
            for (int j=0; j<n; j++) {
                iss >> count;
                M[i][j] = count;
            }
        }
        ifs.close();
        
        partial_ranking.clear();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (M[i][j] < M[j][i])
                    partial_ranking.insert(pair<int, int>(i,j));
            }
        }
        this->partial_rankings.push_back(partial_ranking);
    }
}


RepeatChoice::RepeatChoice (const int &_n, 
	const std::vector<std::vector<std::vector<int> > >& matrice) {
    vector<vector<int> > M;
    set<pair<int, int> > partial_ranking;
    
    this->n = _n;
    
    for (vector<vector<vector<int> > >::const_iterator it = matrice.begin(); 
		it != matrice.end(); it++) {
        M = *it;
        partial_ranking.clear();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (M[i][j] < M[j][i])
                    partial_ranking.insert(pair<int, int>(i,j));
            }
        }
        this->partial_rankings.push_back(partial_ranking);
    }

}


vector<int>
RepeatChoice::process() {
    set<pair<int, int> > ranked_pairs, unranked_pairs, partial_ranking, 
		random_ranking;
    vector<int> ranking_result, random_order, unused_partials;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,this->n);
    pair<int, int> ranked_pair, unranked_pair1, unranked_pair2;
    int random_val;
    
    auto t0 = std::chrono::high_resolution_clock::now();

    //generate a ranking with no object ranked
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++)
            unranked_pairs.insert(pair<int, int>(i, j));
    }
    
    //generate random order
    while (random_order.size() < n) {
        random_val = distribution(generator);
        if (std::find(random_order.begin(), random_order.end(), random_val) 
			== random_order.end())
            random_order.push_back(random_val);
    }
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) {
            random_ranking.insert(pair<int, int>(i, j));
        }
    }
    
    for (int i=0; i<this->partial_rankings.size(); i++)
        unused_partials.push_back(i);
    
    //while unranked object pair is not empty
    while ((!unranked_pairs.empty()) && (!unused_partials.empty())) {
        //find a random partial ranking
        random_val = distribution(generator);
        random_val = random_val % unused_partials.size();
        partial_ranking = this->partial_rankings[unused_partials[random_val]];
    
        //combine it with the current ranking
        for (set<pair<int,int> >::iterator it = partial_ranking.begin(); 
			it != partial_ranking.end(); it ++) {
            unranked_pair1 = *it;
            unranked_pair2 = pair<int, int> (it->second, it->first);
            if ((unranked_pairs.find(unranked_pair1) != unranked_pairs.end()) 
				|| (unranked_pairs.find(unranked_pair2) != unranked_pairs.end())) {
                ranked_pairs.insert(unranked_pair1);
                unranked_pairs.erase(unranked_pair1);
                unranked_pairs.erase(unranked_pair2);
            }
        }
        
        unused_partials.erase(unused_partials.begin()+random_val);
    }
    
    //merge the current ranking with the random order
    for (set<pair<int, int>>::iterator it = random_ranking.begin(); 
		it != random_ranking.end(); it ++) {
        unranked_pair1 = *it;
        unranked_pair2 = pair<int, int> (it->second, it->first);
        if ((unranked_pairs.find(unranked_pair1) != unranked_pairs.end()) 
			|| (unranked_pairs.find(unranked_pair2) != unranked_pairs.end())) {
            ranked_pairs.insert(unranked_pair1);
            unranked_pairs.erase(unranked_pair1);
            unranked_pairs.erase(unranked_pair2);
        }
    }
    
    //organize the ranked pairs
    ranking_result = vector<int> (n);
    for (int i=0; i<n; i++)
        ranking_result[i] = 0;
    
    //count how many elements smaller than i
    for (set<pair<int, int> >::iterator it = ranked_pairs.begin(); 
		it != ranked_pairs.end(); it ++) {
        ranking_result[it->second]++;
    }
    
    vector<set<int> > tmp(n);
    for (int i=0; i<n; i++) {
        tmp[ranking_result[i]].insert(i);
    }
    
    ranking_result.clear();
    
    for (int i=0; i<n; i++) {
        set<int> result = tmp[i];
        for (set<int>::iterator it = result.begin(); it != result.end(); it++)
            ranking_result.push_back(*it);
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
    
    //cout << "RepeatChoice time: " << dt << endl;
    /*
    for (int i=0; i<n; i++)
        cout << ranking_result[i] << " < ";
    cout << endl;
    */
    return ranking_result;
    
}