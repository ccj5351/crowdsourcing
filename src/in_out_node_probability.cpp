#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <map> /*ordered map*/
#include <unordered_map> /*unordered map*/
#include <set>
#include <algorithm> /*Permutation*/
#include <functional>
#include <array> /* array */
#include <cmath> /*pow etc*/
#include <cstdio>/*printf*/

using namespace std;


#define _IN_NODE_OUT_NODE_PROB_
#ifndef _IN_NODE_OUT_NODE_PROB_
int main()
{
	MatrixXd m(2, 2);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;
	m(1, 1) = m(1, 0) + m(0, 1);
	std::cout << m << std::endl;
}
#endif