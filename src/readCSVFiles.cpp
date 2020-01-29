/**
* @file: readCSVFiles.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#include "readCSVFiles.hpp"
using namespace std;
//using namespace 


/*see reference: stackoverflow
 * How can I read and parse CSV files in C++?
 *  http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c#
 */

std::string const &  CSVRow::operator[](std::size_t index) const{
		return v_data[index];
	}

std::size_t CSVRow::size() const{
		return v_data.size();
	}

bool CSVRow::isEmpty(){ return v_data.empty(); }

void CSVRow::readNextRow(std::istream & instream){
		// a line
		std::string line;
		// Extracts characters from "is" and stores them 
		// into "line" until the newline character, '\n'.
		// that is, without explicit delimitation character,
		// "\n" is the default delimitation character therefore.
		std::getline(instream, line);
		std::stringstream   lineStream(line);
		std::string         cell;
		v_data.clear();
		// firstly, we have got a line,
		// then, we extract each word in the line.
		// and those words are separated by comma ",".
		while (std::getline(lineStream, cell, ',')){
			v_data.push_back(cell);
		}
}


/* Stream extraction and insertion:
 * The overloads of operator >> and operator<< that take 
 * a std::istream& or std::ostream& as the left hand argument 
 * are known as insertion and extraction operators.
 * Since they take the user - defined type as the right 
 * argument(b in a@b), they must be implemented as non - members.
 * see http://en.cppreference.com/w/cpp/language/operators
 */
std::istream & operator >> (std::istream & is, CSVRow & data){
	// read data from stream;
	data.readNextRow(is);
	return is;
}

std::vector<CSVRow> readCSVFiles(const string & fname){
	std::ifstream  infile;
	infile.open(fname.c_str());
	std::vector<CSVRow> v_rows;
	if (infile.is_open()){
		while (infile.good()){
			CSVRow row;
			infile >> row;
			if (!row.isEmpty())
				v_rows.push_back(row);
		}
	}
	else
		cout << "Error opening file " << fname << endl;
	infile.close();
	return v_rows;
}

std::vector<vector<double> >
readCSVFile(const string & fname) {
	vector<vector<double> > result;
	vector<string> data;
	int n;

	data = readCSVFiles(fname)[0].v_data;

	n = (int)data.size();

	n = (int)sqrt(n);

	for (int i = 0; i < n; i++) {
		vector<double> row;
		for (int j = 0; j < n; j++)
			row.push_back(stod(data[i*n + j]));
		result.push_back(row);
	}

	return result;
}

bool isCircle(const string & r1, const string& r2, 
	const string & r3, const string & r4){
	// comparing string 1 with string 2
	if ((r1.compare(r2) != 0) || (r2.compare(r3) != 0)
		|| (r3.compare(r4) != 0))
		return false;
	else
		return true;
}

bool isCircle(const char & r1, const char & r2,
	const char & r3, const char & r4){
	// comparing string 1 with string 2
	if ((r1 != r2) || (r2 != r3)
		|| (r3 != r4))
		return false;
	else
		return true;
}