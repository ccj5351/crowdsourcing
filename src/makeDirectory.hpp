/**
* @file: makeDirectory.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef _ReadAndMakeDirectory_H
#define _ReadAndMakeDirectory_H


// boost filesystem
// For new code, defining "BOOST_FILESYSTEM_NO_DEPRECATED" before 
// including filesystem headers is strongly recommended. 
// This prevents inadvertent use of old features, particularly legacy function names, 
// that have been replaced and are going to go away in the future.
#define BOOST_FILESYSTEM_NO_DEPRECATED 
                                      
#include <boost/filesystem/operations.hpp>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
namespace bf = boost::filesystem;

// to create a directory
// returns true if successful;
// otherwise, returns false;
bool MakeDir(const std::string & dir_path);


bool find_file( const bf::path & dir_path,     // in this directory,
                const std::string & file_name, // search for this name,
                bf::path & path_found          // placing path here if found
);

void GetDirList(const std::string & directory, std::vector<std::string> * dirlist);

// to get the names of all the image files in the current path
void GetFileList(const std::string & directory, std::vector<std::string> * filelist);

// overload + 1
void GetFileList(const string& directory, 
	const string & fileextension, // e.g., = ".txt"
	vector<string>* filelist);

// 
/***************************************
*                                      *
* Examples: MakeDir(.), etc            *
*                                      *
***************************************/
/*
 string resultDir = "E:/ComputerVisionCCJ/Epitome/imageDatabase/result";
 MakeDir(resultDir);
 vector<string> categories;
 string databaseDir = "E:/ComputerVisionCCJ/Epitome/imagedata/train";
 GetDirList(databaseDir, &categories);
 vector<string> filelist;
 GetFileList(categories[0], &filelist);
*/

#endif