/**
* @file: main_renameFiles.cpp
* @brief: The main function for the experiments of our paper.
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

/* Usage:
 * What it can do:
    It will extract some characters, between first flag, included, and 
	last flag, excluded, in the filename, and rename it using the 
	extracted characters.

 * Input parameters to the execute file:
 * string baseAddress = argv[1]; // e.g., == "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/GNU-Figures/SA_SATD_2/sa-satd-tt/v10-w20-rt-0.8/";
 * string s_first = argv[2]; // e.g., == "Var";
 * string s_last = argv[3]; // e.g., == "-Smooth";
 * string s_file_extension = argv[4]; // e.g., == ".txt"
 */
#define _RENAME_FILES_BETWEEN_CPP_ /*comment this line, or*/
#ifndef _RENAME_FILES_BETWEEN_CPP_ /*change ifndef to ifdef, to make this main function work.*/
#define _CRT_SECURE_NO_WARNINGS
#include "makeDirectory.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h> /*use rename (or std::rename)*/
using namespace std;


int main(int argc, char * argv[]){
	
	string baseAddress = argv[1]; // e.g., == "E:/OpenCVProjects_CCJ/CrowdSourcing2/ICDE-2017/GNU-Figures/SA_SATD_2/sa-satd-tt/v10-w20-rt-0.8/";
	string s_first = argv[2]; // e.g., == "Var";
	string s_last = argv[3]; // e.g., == "-Smooth";
	string s_file_extension = argv[4]; // e.g., == ".txt"

	vector<string> filelist;
	GetFileList(baseAddress, s_file_extension, &filelist);
	for (int i = 0, s = filelist.size(); i < s; ++i){
		string name1 = filelist[i];
		std::size_t found1 = name1.find(s_first);
		std::size_t found2 = name1.find(s_last);
		if (found1 != std::string::npos && found2 != std::string::npos){
			string name2(name1, found1, found2- found1);
			name2 = name2 + s_file_extension;
			std::rename((baseAddress + name1).c_str(), (baseAddress + name2).c_str());
		}
			
	}
	
	return 0;
}
#endif