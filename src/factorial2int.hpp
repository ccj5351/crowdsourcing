/**
* @file: factorial2int.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_FACTORIAL_TO_INT_H_
#define __HEADER__CROWD_SOURCING_FACTORIAL_TO_INT_H_

typedef unsigned long long ullong;

#include <iostream>
#include <cstdio>
#include <numeric> /*std::accumulate*/
#include <string>
#include <functional> /*std::begin(),std::end() */
#include "MyUtility.hpp"

using namespace std;
class Fact_Int_Mapping {
	/*
	 ********************
	 * What it can do: **
	 ********************
	 * Mapping the combinations of n numbers into integers.
	 * E.g., the combinations of five numbers (1, 2, 3, 4, and 5),
	 * are mapped into integers in the range of [0, 5! - 1] = [0, 119].
	 */

	/*
	 *********************
	 * How it is finished?
	 ********************
	 * For n numbers, the combination is represented
	 * as : (x(0), x(1), x(2), ..., x(n-2), x(n-1)).
	 * For example, 5 numbers: (32415).

	 * The WEIGHT for x(i):
	 * the factorial of the number of digits
	 * on the right of that digit x(i). E.g., the weight of 3 shown
	 * above is 4!, since there are 4 numbers on its right. The weight
	 * of 2 is 3!, the weight of 4 is 2!, and the weight of 1 is 1!.
	 * The last bit is not considered in the mapping, since it will be
	 * uniquely determined when all the other bits are fixed.

	 * The VALUE for each digit:
	 * it depends how many numbers are less than the current digit.
	 * We define the Right-Sub-String of x(i) as the string
	 * containing the elements on the right of x(i), exclusive x(i).
	 * Thus, encode(32415)=2*4!+1*3!+1*2!+0*1!=56.
	 */
private:
	// how many numbers or objects to construct the combination.
	// n_bit <= 20;
	const int n_bit;
	// Since factorial(20) = 2.4329e+18, and 2^64 = 1.8447e+19,
	// thus, if we want to represent the combinations of 20 objects, 
	// we need 64-bit binary integers to finish this 1-to-1 mapping.
	bool isLeq64bits; // Larger than or equal to 64 bits, or 8-byte.

public:
	bool isDisplay;
	Fact_Int_Mapping(const int & n,
		const bool & isDisplay_) :
		n_bit(n), isDisplay(isDisplay_){
		isLeq64bits = (sizeof(ullong) >= 8) ?
			true : false;
		if (n > 20)
			std::cout << "Error! The mapping does not work for "
			<< "the combination of more than 20 objects.\n";
		if (!isLeq64bits)
			std::cout << "Error! There are not 8 bytes.\n";
	};

	// the array permutation = {0, 1, 2, 3, ..., n-1};
	// E.g., (2, 1, 3, 4, 0) ==> 57;  
	ullong fact2int(int * permutation, std::ofstream  & of);

	// permutation = (x0, x1, x2, x3, x4)
	// = (0, 1, 2, 3, 4).
	// mapping 5! combinations to integers within
	// [0, 5! -1] = [0,119].
	// E.g., (2, 1, 3, 4, 0) ==> 57;  
	int fact5_2int(int permutation);


	/* Inverse Mapping:
	 * from integers within
	 * [0, n! -1], backwards n! combinations.
	 * For example, mapping the integer y in the range [0, 5!- 1]
	 * = [0, 119] to the combination C(x0, x1, x2, x3, x4)
	 * of 0, 1, 2, 3, and 4.
	 * For each y, based on the formula:
	 * y=(x0)*4!+(x1)*3!+(x2)*2!+(x3)*1!,
	 * x0 could be 0, 1, 2, 3, or 4, since there are up to
	 * four numbers less than x0 in its Right-Sub-String.
	 * x1 could be 0, 1, 2, 3.
	 * x2 could be 0, 1, 2.
	 * x3 could be 0, 1.
	 * do the following calculation to get:
	 * e.g., y = 57 = 4!*x0+3!*x1+2!*x2+1!*x3+0*x4
	 *         = 24*x0+6*x1+2*x2+1*x3.
	 * x3 = (24*x0+6*x1+2*x2+1*x3) % 2 = y%2 = 57 % 2 =1,
	 * y = (24*x0+6*x1+2*x2+1*x3) /2 = 12*x0+3*x1+1*x2 = y/2 = 28;
	 * x2 = (12*x0+3*x1+1*x2)%3 = y%3 = 28 % 3 =1,
	 * y = (12*x0+3*x1+1*x2)/3 = 4*x0+x1 = y/3 = 9;
	 * x1 = (4*x0+x1)%4 = y%4 = 9 % 4 = 1,
	 * y = (4*x0+x1)%4 = x0 = y/4 = 2;
	 * x0 = y = 2.
	 * Now we get (x0, x1, x2, x3, x4) = (2, 1, 1, 1, x4).
	 * Adjust digit x(i), if its left neighbors x(0 : i-1) is less
	 * then or equal to x(i).
	 * For this example, we have got (x0, x1, x2, x3, x4)=
	 * (2, 1, 1, 1, x4), How to do the adjustment?

	 * Starting from the leftmost digit x(0).
	 * Since x(0) = 2 is already the leftmost digit, thus
	 * it will be left unchanged.
	 *
	 * For x(1)= 1, since its left neighbors x(0), is larger than x(1),
	 * which will not affect how many digits are less than x(1) in the
	 * Right-Sub-String of x(1). Thus leave x(1) unchanged.
	 *
	 * For x(2), since x(1) <= x(2), thus x(2) += 1 = 2; the new value
	 * of x(2) >= x(0), thus x(2) += 1 = 3;
	 *
	 * For x(3) = 1, since the old values x(2) = 1, x(2) <= x(3), thus
	 * x(3) += 1 = 2. The new x(3) >= x(1)(the old value) = 1, thus
	 * x(3) += 1 = 3. The new x(3) >= x(0)(the old value) = 2, thus
	 * x(3) += 1 = 4.
	 *
	 * For the last digit x(4), since the sum of those numbers are fixed,
	 * therefore, x(4) = sum(0: 4) - sum(x(0), x(1), x(2), x(3)) = 0.
	 * Finally, we have the inverse mapping:
	 * 57 ==> (2, 1, 3, 4, 0).
	 */
	// E.g., 57 ==> (2, 1, 3, 4, 0)
	int * int2fact(const ullong & gene,
		std::string combi);

	vector<int>  int2fact2(const ullong & gene,
		std::string combi, std::ofstream  & of);

	// inverse mapping:
	// that is : mapping integers within
	// [0, 5! -1] = [0,119], backwards
	// 5! combinations.
	// E.g., 57 ==> (2, 1, 3, 4, 0)  
	int int2fact5(const int & gene);
	
};

#endif
