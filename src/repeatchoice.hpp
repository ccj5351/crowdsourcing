/**
* @file: repeatchoice.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef repeatchoice_hpp
#define repeatchoice_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <utility>
#include <set>

class RepeatChoice {
public:
    RepeatChoice() {};
    RepeatChoice(const int &_n, const std::vector<std::string> &_files);
    RepeatChoice(const int &_n, const std::vector<std::vector<std::vector<int> > >& matrice);
    std::vector<int> process ();    //generate a total ranking
    
    int n;  //# of objects
    std::vector<std::set<std::pair<int, int> > > partial_rankings;  //each partial ranking is a set of pairs, i < j
    std::vector<std::string> files; //the files that contain the partial ranking
};

#endif /* repeatchoice_hpp */
