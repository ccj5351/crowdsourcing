/**
* @file: readCSVFiles.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_READ_CSV_FILE_H_
#define __HEADER__CROWD_SOURCING_READ_CSV_FILE_H_
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;
//using namespace 
class CSVRow{

public:
	std::vector<std::string>    v_data;

public:
	//deconstructor
	~CSVRow(){ vector<string>().swap(v_data); }
	std::string const & operator[](std::size_t index) const;
	std::size_t size() const;
	void readNextRow(std::istream & instream);
	bool isEmpty();
};

/* Stream extraction and insertion:
* The overloads of operator >> and operator<< that take
* a std::istream& or std::ostream& as the left hand argument
* are known as insertion and extraction operators.
* Since they take the user - defined type as the right
* argument(b in a@b), they must be implemented as non - members.
* see http://en.cppreference.com/w/cpp/language/operators
*/
std::istream & operator >> (std::istream & is, CSVRow & data);

std::vector<CSVRow> readCSVFiles(const string & fname);
std::vector<vector<double>> readCSVFile(const string & fname);
bool isCircle(const char & r1, const char & r2,
	const char & r3, const char & r4);
bool isCircle(const string & r1, const string& r2,
	const string & r3, const string & r4);

#endif