/**
* @file: ktaub.hpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

// see https://afni.nimh.nih.gov/pub/dist/src/ktaub.c
/*-------------------------------------------------------------------------*/
/* This file contains code to calculate Kendall's Tau in O(N log N) time in
* a manner similar to the following reference:
*
* A Computer Method for Calculating Kendall's Tau with Ungrouped Data
* William R. Knight Journal of the American Statistical Association,
* Vol. 61, No. 314, Part 1 (Jun., 1966), pp. 436-439
*
* Copyright 2010 David Simcha
*
* License:
* Boost Software License - Version 1.0 - August 17th, 2003
*
* Permission is hereby granted, free of charge, to any person or organization
* obtaining a copy of the software and accompanying documentation covered by
* this license (the "Software") to use, reproduce, display, distribute,
* execute, and transmit the Software, and to prepare derivative works of the
* Software, and to permit third-parties to whom the Software is furnished to
* do so, all subject to the following:
*
* The copyright notices in the Software and this entire statement, including
* the above license grant, this restriction and the following disclaimer,
* must be included in all copies of the Software, in whole or in part, and
* all derivative works of the Software, unless such copies or derivative
* works are solely in the form of machine-executable object code generated by
* a source language processor.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*/
/*-------------------------------------------------------------------------*/
/**---------------- Adapted for AFNI by RWCox - April 2010 ---------------**/
/*-------------------------------------------------------------------------*/
#ifndef __HEADER__CROWD_SOURCING_KENDALL_TAU_DISTANCE_H_
#define __HEADER__CROWD_SOURCING_KENDALL_TAU_DISTANCE_H_
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdio>
#include <assert.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>      // std::fstream

#ifdef SOLARIS_OLD      /* 'if' might be overkill, but just to be minimal... */
/*                               10 May 2010 [rickr] */
#include "machdep.h"
#else
#include <stdint.h>
#endif

/*------------------------------------------------------*/
/* Sorts in place, returns the bubble sort distance between
* the input array and the sorted array.
*/
static int insertionSort(float *arr, int len);

/*------------------------------------------------------*/

static int merge(float *from, float *to, int middle, int len);

/*---------------------------------------------------*/
/* Sorts in place, returns the bubble sort distance
* between the input array and the sorted array.
*/

static int mergeSort(float *x, float *buf, int len);

/*-------------------------------------------------------*/

static int getMs(float *data, int len);

/*-------------------------------------------------------*/
/* This function calculates the Kendall correlation tau_b.
* The arrays arr1 should be sorted before this call,
* and arr2 should be re-ordered in lockstep.  This can be
* done by calling qsort_floatfloat(len,arr1,arr2) for example.
* Note also that arr1 and arr2 will be modified, so if they
* need to be preserved, do so before calling this function.
*/

using namespace std;
float kendallNlogN(float *arr1, float *arr2, int len);

/*----------------------------------------------------------*/
/* This function uses a simple O(N^2) implementation.
* It probably has a smaller constant and therefore
* is useful in the small N case, and is also
* useful for testing the relatively complex O(N log N)
* implementation.
*/

float  kendallSmallN(float *arr1, float *arr2, int len);

double kendallSmallN(
	const int * arr1,
	const string & a1ID,
	const int *arr2,
	const string & a2ID,
	const int &  len,
	const bool & isDisplay);

enum CompareResult : int { GREATER = 1, LESS = 2, EQUAL = 3 };
inline CompareResult compareRank(const int & r1,
	const int & r2, const bool & isDisplay){
	if (r1 > r2){
		return GREATER;
		if (isDisplay)
			std::cout << r1 << " > " << r2 << std::endl;
	}
	else if (r1 < r2){
		return LESS;
		if (isDisplay)
			std::cout << r1 << " < " << r2 << std::endl;
	}
	else{
		return EQUAL;
		if (isDisplay)
			std::cout << r1 << " = " << r2 << std::endl;
	}
}

inline void printCompareResult(const int & r1,
	const int & r2,
	const CompareResult & re){
	if (re == GREATER){
		std::cout << r1 << " > " << r2 << ", ";
	}
	else if (re == LESS){
		std::cout << r1 << " < " << r2 << ", ";
	}
	else{
		std::cout << r1 << " = " << r2 << ", ";
	}
}

// A complete testing of Kendall Tau distance for 
// different ranking lists.
int execute_Ktau_main2();

// only for Kendall Tau distance for small length ranking lists.
void execute_Ktau_main1(const bool & isDisplay,
	std::ofstream  & of);

#endif