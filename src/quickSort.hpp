/**
* @file: quickSort.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#ifndef __HEADER__CROWD_SOURCING_QUICKSORT_H_
#define __HEADER__CROWD_SOURCING_QUICKSORT_H_

#include <iostream>
#include <algorithm>
//see QuickSort code: http://ideone.com/ImMUH4
void print(int *a, int n);

int partition(int *arr, const int left, const int right);

void quicksort(int *arr, const int left, const int right, const int sz);
#endif