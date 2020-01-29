/**
* @file: transition.hpp
* @brief: to compute the transitivity closure.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef transition_hpp__
#define transition_hpp__

#include <stdio.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

//read matrix from the file, and generate the transition matrix.
std::vector<std::vector<double> > 
generateTransitionMatrix (const std::string &path, const int &n);    


std::vector<std::vector<double> > 
generateTransitionMatrix (const std::vector<std::vector<double> >& matrix);

std::vector<std::vector<double> >
generateTransitionMatrix_Eigen(const std::vector<std::vector<double>>& matrix);

#endif /* transition_hpp */
