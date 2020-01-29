/**
* @file: factorial2int.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#include "factorial2int.hpp"

// the array permutation = {0, 1, 2, 3, ..., n-1};
// E.g., (2, 1, 3, 4, 0) ==> 57;  
ullong Fact_Int_Mapping::fact2int(int * permutation, std::ofstream  & of){
	// input: 
	// combination of n objects array: x[]:
	// shown as (x(0), x(1), x(2), ..., x(n-2), x(n-1)).
	int * x = new int[n_bit];
	ullong ret = 0;
	// initialization.
	for (int i = 0; i < n_bit; ++i){
		x[i] = permutation[i];
	}

	std::string combi =
		std::accumulate(x, x + n_bit, std::string{},
		[](const std::string & a, int b) {
		return a.empty() ? std::to_string(b)
			: a + '-' + std::to_string(b);
	});

	for (int i = 0; i < n_bit - 1; // the last bit is not involved.
		++i){
		// In the Right-Sub-String of the current bit x(i):
		// count how many bits are less than it.
		for (int j = i + 1; j < n_bit - 1; ++j){
			x[j] = x[j] > x[i] ? x[j] - 1 : x[j];
		}
		ret += factorial(n_bit - i - 1) * (x[i]);
	}
	if (isDisplay){
		std::cout << "Combination " << combi
		<< " is mapped to integer " << ret << ".\n";
		of << "Combination " << combi
			<< " is mapped to integer " << ret << ".\n";
	}
	delete[] x;
	return ret;
}

// permutation = (x0, x1, x2, x3, x4)
// = (0, 1, 2, 3, 4).
// mapping 5! combinations to integers within
// [0, 5! -1] = [0,119].
// E.g., (2, 1, 3, 4, 0) ==> 57;  
int Fact_Int_Mapping::fact5_2int(int permutation){
	int x0 = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0;
	int ret = 0;
	// initialize
	x0 = permutation / 10000;
	permutation -= 10000 * x0;
	x1 = permutation / 1000;
	permutation -= 1000 * x1;
	x2 = permutation / 100;
	permutation -= 100 * x2;
	x3 = permutation / 10;
	permutation -= 10 * x3;
	x4 = permutation;
	string combi = to_string(x0) + "-" + to_string(x1) + "-"
		+ to_string(x2) + "-" + to_string(x3) + "-" + to_string(x4);
	// for the leftmost bit x0.
	// its Right-Sub-String is (x1, x2, x3, x4)
	// but the last bit x4 is not involved.
	ret += 24 * (x0);
	x1 = x1 > x0 ? x1 - 1 : x1;
	x2 = x2 > x0 ? x2 - 1 : x2;
	x3 = x3 > x0 ? x3 - 1 : x3;
	ret += 6 * (x1);
	// for the bit x1.
	// its Right-Sub-String is (x2, x3, x4).
	// but the last bit x4 is not involved.
	x2 = x2 > x1 ? x2 - 1 : x2;
	x3 = x3 > x1 ? x3 - 1 : x3;
	ret += 2 * (x2);
	// for the bit x2.
	// its Right-Sub-String is (x3, x4).
	// but the last bit x4 is not involved.
	x3 = x3 > x2 ? x3 - 1 : x3;
	ret += x3;

	if (isDisplay){
		std::cout << "Combination " << combi
			<< " is mapped to integer " << ret << ".\n";
	}
	// done! return the result.
	return ret;
}


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
int * Fact_Int_Mapping::int2fact(const ullong & gene,
	std::string combi){
	// Output:
	// combination of n objects array: x[]:
	// shown as (x(0), x(1), x(2), ..., x(n-2), x(n-1)).
	int * x = new int[n_bit];
	// array b is used to count how many digits are 
	// less than or equal to the current digit.
	int * b = new int[n_bit];
	// initialization.
	for (int i = 0; i < n_bit; ++i){
		x[i] = 0;
		b[i] = 0;
	}
	ullong temp_gene = gene;

	/*Calculate the position in Right-Sub_string*/
	/*exclude the first element x(0), and the last
	* element x(n-1).
	*/

	/* E.g., y = 57 = 4!*x0 + 3!*x1 + 2!*x2 + 1!*x3 + 0 * x4
	*         = 24 * x0 + 6 * x1 + 2 * x2 + 1 * x3.
	* x3 = (24 * x0 + 6 * x1 + 2 * x2 + 1 * x3) % 2
	*    = y % 2 = 57 % 2 = 1,
	* y  = (24 * x0 + 6 * x1 + 2 * x2 + 1 * x3) / 2
	*    = 12 * x0 + 3 * x1 + 1 * x2 = y / 2 = 28;
	* x2 = (12 * x0 + 3 * x1 + 1 * x2) % 3 = y % 3 = 28 % 3 = 1,
	* y  = (12 * x0 + 3 * x1 + 1 * x2) / 3
	*    = 4 * x0 + x1 = y / 3 = 9;
	* x1 = (4 * x0 + x1) % 4 = y % 4 = 9 % 4 = 1,
	* y  = (4 * x0 + x1) % 4 = x0 = y / 4 = 2;
	* x0 = y = 2.
	*/
	for (int i = n_bit - 2; i > 0; --i){
		// Ordering is :
		// starting from the the last second element, 
		// to the first second element.
		// that is, from x(n-2) to x(1). 
		// since the last element x(n-1), and the first x(0)
		// are not involved in the inverse mapping.
		x[i] = temp_gene % (n_bit - i);
		temp_gene /= (n_bit - i);
	}

	// for the first element
	x[0] = temp_gene;

	/* Currently, we have obtained the value of each digit,
	* based on their Right-Sub-Strings.
	* We further need to adjust their values,
	* by considering the global string. Specifically,
	* for each digit x(i), if its left neighbors x(0 : i-1) are less
	* than or equal to x(i), adjust x(i) by adding 1, i.e.,
	* x(i) ++, respectively.
	*/
	// the outer for-loop, staring from the first second element,
	// ending up with the last second element.
	for (int i = 1; i < n_bit - 1; ++i){
		// the inner for-loop, counting how many elements
		// (i.e., digits) are less than the current element, by 
		// searching through its left neighbors.
		for (int j = i - 1; j >= 0; --j){
			b[i] += (x[i] + b[i] >= x[j]) ? 1 : 0;
		}
	}

	for (int i = 1; i < n_bit - 1; ++i){
		x[i] += b[i];
	}

	// for the last element x(n-1),
	// which is not involved in the mapping.
	ullong sum = 0, sum1 = 0;
	for (int i = 0; i < n_bit; ++i){
		sum += i;
	}
	for (int i = 0; i < n_bit - 1; ++i){
		sum1 += x[i];
	}

	x[n_bit - 1] = sum - sum1;
	//see http://en.cppreference.com/w/cpp/algorithm/accumulate;
	combi = std::accumulate(x, x + n_bit, std::string{},
		[](const std::string & a, int b) {
		return a.empty() ? std::to_string(b)
			: a + '-' + std::to_string(b);
	});
	if (isDisplay){
		std::cout << "Integer " << gene
		<< " is mapped to combination " << combi << ".\n";
	}
	return x;
}


vector<int>  Fact_Int_Mapping::int2fact2(const ullong & gene,
	std::string combi, 
	std::ofstream  & of){
	// Output:
	// combination of n objects array: x[]:
	// shown as (x(0), x(1), x(2), ..., x(n-2), x(n-1)).
	vector<int> x(n_bit);
	// array b is used to count how many digits are 
	// less than or equal to the current digit.
	int * b = new int[n_bit];
	// initialization.
	for (int i = 0; i < n_bit; ++i){
		x[i] = 0;
		b[i] = 0;
	}
	ullong temp_gene = gene;

	/*Calculate the position in Right-Sub_string*/
	/*exclude the first element x(0), and the last
	* element x(n-1).
	*/

	/* E.g., y = 57 = 4!*x0 + 3!*x1 + 2!*x2 + 1!*x3 + 0 * x4
	*         = 24 * x0 + 6 * x1 + 2 * x2 + 1 * x3.
	* x3 = (24 * x0 + 6 * x1 + 2 * x2 + 1 * x3) % 2
	*    = y % 2 = 57 % 2 = 1,
	* y  = (24 * x0 + 6 * x1 + 2 * x2 + 1 * x3) / 2
	*    = 12 * x0 + 3 * x1 + 1 * x2 = y / 2 = 28;
	* x2 = (12 * x0 + 3 * x1 + 1 * x2) % 3 = y % 3 = 28 % 3 = 1,
	* y  = (12 * x0 + 3 * x1 + 1 * x2) / 3
	*    = 4 * x0 + x1 = y / 3 = 9;
	* x1 = (4 * x0 + x1) % 4 = y % 4 = 9 % 4 = 1,
	* y  = (4 * x0 + x1) % 4 = x0 = y / 4 = 2;
	* x0 = y = 2.
	*/
	for (int i = n_bit - 2; i > 0; --i){
		// Ordering is :
		// starting from the the last second element, 
		// to the first second element.
		// that is, from x(n-2) to x(1). 
		// since the last element x(n-1), and the first x(0)
		// are not involved in the inverse mapping.
		x[i] = temp_gene % (n_bit - i);
		temp_gene /= (n_bit - i);
	}

	// for the first element
	x[0] = temp_gene;

	/* Currently, we have obtained the value of each digit,
	* based on their Right-Sub-Strings.
	* We further need to adjust their values,
	* by considering the global string. Specifically,
	* for each digit x(i), if its left neighbors x(0 : i-1) are less
	* than or equal to x(i), adjust x(i) by adding 1, i.e.,
	* x(i) ++, respectively.
	*/
	// the outer for-loop, staring from the first second element,
	// ending up with the last second element.
	for (int i = 1; i < n_bit - 1; ++i){
		// the inner for-loop, counting how many elements
		// (i.e., digits) are less than the current element, by 
		// searching through its left neighbors.
		for (int j = i - 1; j >= 0; --j){
			b[i] += (x[i] + b[i] >= x[j]) ? 1 : 0;
		}
	}

	for (int i = 1; i < n_bit - 1; ++i){
		x[i] += b[i];
	}

	// for the last element x(n-1),
	// which is not involved in the mapping.
	ullong sum = 0, sum1 = 0;
	for (int i = 0; i < n_bit; ++i){
		sum += i;
	}
	for (int i = 0; i < n_bit - 1; ++i){
		sum1 += x[i];
	}

	x[n_bit - 1] = sum - sum1;
	//see http://en.cppreference.com/w/cpp/algorithm/accumulate;
	combi = std::accumulate(x.begin(), x.end(), std::string{},
		[](const std::string & a, int b) {
		return a.empty() ? std::to_string(b)
			: a + '-' + std::to_string(b);
	});
	if (isDisplay){
		std::cout << "Integer " << gene
		<< " is mapped to combination " << combi << ".\n";
		of << "Integer " << gene
			<< " is mapped to combination " << combi << ".\n";
	}
	delete [] b;
	return x;
}


// inverse mapping:
// that is : mapping integers within
// [0, 5! -1] = [0,119], backwards
// 5! combinations.
// E.g., 57 ==> (2, 1, 3, 4, 0)  
int Fact_Int_Mapping::int2fact5(const int & gene){
	int ret = 0;
	int x0 = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0;
	int temp_gene = gene;
	/* gene = 4!*x0 + 3!*x1 + 2!*x2 + 1!*x3 + 0 * x4*/
	x3 = temp_gene % 2;
	temp_gene /= 2;

	/*gene=12*x0+3*x1+x2*/
	x2 = temp_gene % 3;
	temp_gene /= 3;

	/*gene=4*x0+x1;*/
	x1 = temp_gene % 4;
	temp_gene /= 4;

	/*gene=x1;*/
	x0 = temp_gene;
	/*
	if the permutation involves (0, 1, 2, 3, 4):
	do not do the following "++" operation.
	x0++; x1++; x2++; x3++;
	if the permutation involves (1, 2, 3, 4, 5):
	do the following "++" operation.
	x0++; x1++; x2++; x3++;
	*/

	int b1 = 0, b2 = 0, b3 = 0;
	b1 += x1 + b1 >= x0 ? 1 : 0;

	b2 += x2 + b2 >= x1 ? 1 : 0;
	b2 += x2 + b2 >= x0 ? 1 : 0;

	b3 += x3 + b3 >= x2 ? 1 : 0;
	b3 += x3 + b3 >= x1 ? 1 : 0;
	b3 += x3 + b3 >= x0 ? 1 : 0;

	x1 += b1;
	x2 += b2;
	x3 += b3;

	// for the last element x(4),
	// which is not involved in the mapping.
	x4 = (1 + 2 + 3 + 4) - x0 - x1 - x2 - x3;
	ret = x0 * 10000 + x1 * 1000 + x2 * 100 + x3 * 10 + x4;
	string combi = to_string(x0) + "-" + to_string(x1) + "-"
		+ to_string(x2) + "-" + to_string(x3) + "-" + to_string(x4);
	if (isDisplay){
		std::cout << "Integer " << gene
			<< " is mapped to combination " << combi << ".\n";
	}
	return ret;
}