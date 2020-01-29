/**
* @file: ktaub.cpp
* @brief:
* @author: Changjiang Cai, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/
#include "ktaub.hpp"
using namespace std;
/*------------------------------------------------------*/
/* Sorts in place, returns the bubble sort distance between
* the input array and the sorted array.
*/
static int insertionSort(float *arr, int len)
{
	int maxJ, i, j, swapCount = 0;
	/* printf("enter insertionSort len=%d\n",len) ; */

	if (len < 2) { return 0; }

	maxJ = len - 1;
	for (i = len - 2; i >= 0; --i) {
		float  val = arr[i];
		for (j = i; j < maxJ && arr[j + 1] < val; ++j) {
			arr[j] = arr[j + 1];
		}

		arr[j] = val;
		swapCount += (j - i);
	}

	return swapCount;
}

/*------------------------------------------------------*/

static int merge(float *from, float *to, int middle, int len)
{
	int bufIndex, leftLen, rightLen, swaps;
	float *left, *right;

	/* printf("enter merge\n") ; */

	bufIndex = 0;
	swaps = 0;

	left = from;
	right = from + middle;
	rightLen = len - middle;
	leftLen = middle;

	while (leftLen && rightLen) {
		if (right[0] < left[0]) {
			to[bufIndex] = right[0];
			swaps += leftLen;
			rightLen--;
			right++;
		}
		else {
			to[bufIndex] = left[0];
			leftLen--;
			left++;
		}
		bufIndex++;
	}

	if (leftLen) {
#pragma omp critical (MEMCPY)
		memcpy(to + bufIndex, left, leftLen * sizeof(float));
	}
	else if (rightLen) {
#pragma omp critical (MEMCPY)
		memcpy(to + bufIndex, right, rightLen * sizeof(float));
	}

	return swaps;
}

/*---------------------------------------------------*/
/* Sorts in place, returns the bubble sort distance
* between the input array and the sorted array.
*/

static int mergeSort(float *x, float *buf, int len)
{
	int swaps, half;

	/* printf("enter mergeSort\n") ; */

	if (len < 10) {
		return insertionSort(x, len);
	}

	swaps = 0;

	if (len < 2) { return 0; }

	half = len / 2;
	swaps += mergeSort(x, buf, half);
	swaps += mergeSort(x + half, buf + half, len - half);
	swaps += merge(x, buf, half, len);

#pragma omp critical (MEMCPY)
	memcpy(x, buf, len * sizeof(float));
	return swaps;
}

/*-------------------------------------------------------*/

static int getMs(float *data, int len)  /* Assumes data is sorted */
{
	int Ms = 0, tieCount = 0, i;

	/* printf("enter getMs\n") ; */

	for (i = 1; i < len; i++) {
		if (data[i] == data[i - 1]) {
			tieCount++;
		}
		else if (tieCount) {
			Ms += (tieCount * (tieCount + 1)) / 2;
			tieCount = 0;
		}
	}
	if (tieCount) {
		Ms += (tieCount * (tieCount + 1)) / 2;
	}
	return Ms;
}

/*-------------------------------------------------------*/
/* This function calculates the Kendall correlation tau_b.
* The arrays arr1 should be sorted before this call,
* and arr2 should be re-ordered in lockstep.  This can be
* done by calling qsort_floatfloat(len,arr1,arr2) for example.
* Note also that arr1 and arr2 will be modified, so if they
* need to be preserved, do so before calling this function.
*/

float kendallNlogN(float *arr1, float *arr2, int len)
{
	int m1 = 0, m2 = 0, tieCount, swapCount, nPair, s, i;
	float cor;

	/* printf("enter kendallNlogN\n") ; */

	if (len < 2) return (float)0;

	nPair = len * (len - 1) / 2;
	s = nPair;

	tieCount = 0;
	for (i = 1; i < len; i++) {
		if (arr1[i - 1] == arr1[i]) {
			tieCount++;
		}
		else if (tieCount > 0) {
			insertionSort(arr2 + i - tieCount - 1, tieCount + 1);
			m1 += tieCount * (tieCount + 1) / 2;
			s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
			tieCount = 0;
		}
	}
	if (tieCount > 0) {
		insertionSort(arr2 + i - tieCount - 1, tieCount + 1);
		m1 += tieCount * (tieCount + 1) / 2;
		s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
	}

	swapCount = mergeSort(arr2, arr1, len);

	m2 = getMs(arr2, len);
	s -= (m1 + m2) + 2 * swapCount;

	if (m1 < nPair && m2 < nPair)
		cor =
		s / (sqrtf((float)(nPair - m1)) * sqrtf((float)(nPair - m2)));
	else
		cor = 0.0f;

	return cor;
}



// two input arrays arr1 and arr2 should contain the same elements,
// but with different rankings.
double kendallSmallN(const int *arr1,
	const string & a1ID,
	const int *arr2,
	const string & a2ID,
	const int & len,
	const bool &isDisplay)
{
	int s = 0, nPair = len * (len - 1) / 2;
	std::vector<int> v_orderd(arr1, arr1 + len);
	// Sorts the elements in the range[first, last) 
	// into ascending order.
	std::sort(v_orderd.begin(), v_orderd.end());
	/*
	// print out content:
	std::cout << "The ascending sorted array contains:";
	for (std::vector<int>::iterator it = v_orderd.begin();
	it != v_orderd.end(); ++it)
	std::cout << ' ' << *it;
	std::cout << '\n';
	*/
	std::map<int, int> m_rankArry1, m_rankArry2;
	for (int i = 0; i < len; i++) {
		m_rankArry1[arr1[i]] = i;
		m_rankArry2[arr2[i]] = i;
	}
	if (isDisplay){
		std::cout << "Compare " << a1ID << " and "
			<< a2ID << ":\n";
	}

	for (int i = 0; i < len - 1; i++) {
		for (int j = i + 1; j < len; j++) {
			CompareResult t1 = compareRank(m_rankArry1[v_orderd[i]],
				m_rankArry1[v_orderd[j]], false),
				t2 = compareRank(m_rankArry2[v_orderd[i]],
				m_rankArry2[v_orderd[j]], false);

			if (t1 != t2){
				s++;
				if (isDisplay){
					std::cout << "Disagreement " << s << " : ";
					std::cout << a1ID << " gives ";
					printCompareResult(v_orderd[i],
						v_orderd[j], t1);
					std::cout << ", and " << a2ID << " gives ";
					printCompareResult(v_orderd[i],
						v_orderd[j], t2);
					cout << endl;
				}
			}
		}
	}
	double r = double(s) / double(nPair);
	if (isDisplay)
		std::cout << "Kendall Tau Distance is " << s <<
		"/" << nPair << " = " << r << std::endl;
	return r;
}

/*----------------------------------------------------------*/
/* This function uses a simple O(N^2) implementation.
* It probably has a smaller constant and therefore
* is useful in the small N case, and is also
* useful for testing the relatively complex O(N log N)
* implementation.
*/
float kendallSmallN(float *arr1, float *arr2, int len)
{
	int m1 = 0, m2 = 0, s = 0, nPair, i, j;
	float cor;

	/* printf("enter kendallSmallN\n") ; */

	for (i = 0; i < len; i++) {
		for (j = i + 1; j < len; j++) {
			if (arr2[i] > arr2[j]) {
				if (arr1[i] > arr1[j]) {
					s++;
				}
				else if (arr1[i] < arr1[j]) {
					s--;
				}
				else {
					m1++;
				}
			}
			else if (arr2[i] < arr2[j]) {
				if (arr1[i] > arr1[j]) {
					s--;
				}
				else if (arr1[i] < arr1[j]) {
					s++;
				}
				else {
					m1++;
				}
			}
			else {
				m2++;

				if (arr1[i] == arr1[j]) {
					m1++;
				}
			}
		}
	}

	nPair = len * (len - 1) / 2;

	if (m1 < nPair && m2 < nPair)
		cor =
		s / (sqrtf((float)(nPair - m1)) * sqrtf((float)(nPair - m2)));
	else
		cor = 0.0f;

	return cor;
}

// A complete testing of Kendall Tau distance for 
// different ranking lists.
int execute_Ktau_main2(){
	float  a[100], b[100];
	float  smallNCor, smallNCov, largeNCor, largeNCov;
	int i;
	int N;
#define NBOT  10
#define NTOP  2000
#define NSTEP 10

	/* Test the small N version against a few values
	* obtained from the old version in R.
	* Only exercising the small N version because the large
	* N version requires the first array to be sorted
	* and the second to be reordered in lockstep before
	* it's called.
	*/
	{
		float   a1[] = { 1, 2, 3, 4, 5 };
		float   a2[] = { 3, 4, 1, 2, 5 };
		float   b1[] = { 8, 6, 7, 5, 3, 0, 9 };
		float   b2[] = { 3, 1, 4, 1, 5, 9, 2 };
		float   c1[] = { 1, 1, 1, 2, 3, 3, 4, 4 };
		float   c2[] = { 1, 2, 1, 3, 3, 5, 5, 5 };

		printf("kSmall 5 = %g\n", kendallSmallN(a1, a2, 5));
		printf("kSmall 7 = %g\n", kendallSmallN(b1, b2, 7));
		printf("kSmall 8 = %g\n", kendallSmallN(c1, c2, 8));
	}

	/* Now that we're confident that the simple, small N version works,
	* extensively test it against the much more complex and bug-prone
	* O(N log N) version.
	*/
	for (i = 0; i < 100; i++) {
		int j, len;
		for (j = 0; j < 100; j++) {
			a[j] = rand() % 30;
			b[j] = rand() % 30;
		}

		len = rand() % 50 + 50;

		/* The large N version assumes that the first array is sorted.
		* This will usually be made true in R before passing the arrays
		* to the C functions.
		*/
		insertionSort(a, len);

		smallNCor = kendallSmallN(a, b, len);
		largeNCor = kendallNlogN(a, b, len);
		printf("Cor: kSmall = %g  kNlogN = %g\n", smallNCor, largeNCor);
	}

	printf("-- short tests done --\n");

	/* Speed test.  Compare the O(N^2) version,
	* which is very similar to
	* R's current impl, to my O(N log N) version.
	*/
	{
		float  *foo, *bar, *buf;
		int i;
		float  startTime, stopTime;

		foo = (float *)malloc(NTOP * sizeof(float));
		bar = (float *)malloc(NTOP * sizeof(float));
		buf = (float *)malloc(NTOP * sizeof(float));

		for (N = NBOT; N <= NTOP; N += NSTEP){
			for (i = 0; i < N; i++) {
				foo[i] = rand();
				bar[i] = rand();
			}

			startTime = clock();
			smallNCor = kendallSmallN(foo, bar, N);
			stopTime = clock();
			printf("N=%d: slow = %f ms  val = %g\n",
				N, stopTime - startTime, smallNCor);

			startTime = clock();

			/* Only sorting first array.
			* Normally the second one would be
			* reordered in lockstep.
			*/

			mergeSort(foo, buf, N);
			largeNCor = kendallNlogN(foo, bar, N);
			stopTime = clock();
			printf("N=%d: fast = %f ms  val = %g\n",
				N, stopTime - startTime, largeNCor);
		}
	}

	return 0;
}

/*
vector<int>
pE20W5N = { 0, 2, 8, 5, 6, 3, 4, 7, 1, 9, }, // "E20W5";
pE25W5N = { 0, 6, 3, 2, 8, 5, 1, 9, 4, 7, },
pE25W5O = { 0, 2, 3, 6, 5, 8, 1, 4, 9, 7, },
pE30W5N = { 0, 2, 3, 6, 8, 5, 1, 4, 9, 7, },
pE30W5O = { 0, 2, 6, 8, 3, 5, 4, 1, 9, 7, },
pE40W5N = { 0, 2, 3, 6, 5, 8, 1, 4, 9, 7, },
pE40W5O = { 0, 2, 3, 6, 8, 5, 1, 4, 9, 7, },
pE20W7O = { 0, 2, 6, 3, 8, 4, 5, 7, 1, 9, },
pE20W10N = { 0, 2, 8, 3, 6, 5, 9, 4, 7, 1, },
pE20W10O = { 0, 2, 6, 5, 8, 3, 4, 7, 1, 9, },
pAlgo = { 0, 1, 5, 6, 4, 2, 7, 3, 9, 8, };

std::pair<string, vector<int>>
spE20W5N = make_pair("E20W5N", pE20W5N),
spE25W5N = make_pair("E25W5N", pE25W5N),
spE25W5O = make_pair("E25W5O", pE25W5O),
spE30W5N = make_pair("E30W5N", pE30W5N),
spE30W5O = make_pair("E30W5O", pE30W5O),
spE40W5N = make_pair("E40W5N", pE40W5N),
spE40W5O = make_pair("E40W5O", pE40W5O),
spE20W7O = make_pair("E20W7O", pE20W7O),
spE20W10N = make_pair("E20W10N", pE20W10N),
spE20W10O = make_pair("E20W10O", pE20W10O),
spAlgo = make_pair("Algori", pAlgo);
vector<std::pair<string, vector<int>>> v_m_temp = {
spE20W5N,
spE25W5N,
spE25W5O,
spE30W5N,
spE30W5O,
spE40W5N,
spE40W5O,
spE20W7O,
spE20W10N,
spE20W10O,
spAlgo,
};
*/
// only for Kendall Tau distance for small length ranking lists.
void execute_Ktau_main1(const bool & isDisplay,
	std::ofstream  & of){
	int len = 10;
	vector<int>
		pE20W5N = { 0, 2, 8, 5, 6, 3, 4, 7, 1, 9, }, // "E20W5";
		pE25W5N = { 0, 6, 3, 2, 8, 5, 1, 9, 4, 7, },
		pE25W5O = { 0, 2, 3, 6, 5, 8, 1, 4, 9, 7, },
		pE30W5N = { 0, 2, 3, 6, 8, 5, 1, 4, 9, 7, },
		pE30W5O = { 0, 2, 6, 8, 3, 5, 4, 1, 9, 7, },
		pE40W5N = { 0, 2, 3, 6, 5, 8, 1, 4, 9, 7, },
		pE40W5O = { 0, 2, 3, 6, 8, 5, 1, 4, 9, 7, },
		pE20W7O = { 0, 2, 6, 3, 8, 4, 5, 7, 1, 9, },
		pE20W10N = { 0, 2, 8, 3, 6, 5, 9, 4, 7, 1, },
		pE20W10O = { 0, 2, 6, 5, 8, 3, 4, 7, 1, 9, },
		pAlgo = { 0, 1, 5, 6, 4, 2, 7, 3, 9, 8, };

	std::pair<string, vector<int>>
		spE20W5N = make_pair("E20W5N", pE20W5N),
		spE25W5N = make_pair("E25W5N", pE25W5N),
		spE25W5O = make_pair("E25W5O", pE25W5O),
		spE30W5N = make_pair("E30W5N", pE30W5N),
		spE30W5O = make_pair("E30W5O", pE30W5O),
		spE40W5N = make_pair("E40W5N", pE40W5N),
		spE40W5O = make_pair("E40W5O", pE40W5O),
		spE20W7O = make_pair("E20W7O", pE20W7O),
		spE20W10N = make_pair("E20W10N", pE20W10N),
		spE20W10O = make_pair("E20W10O", pE20W10O),
		spAlgo = make_pair("Algori", pAlgo);

	vector<std::pair<string, vector<int>>> v_m_temp = {
		//spE20W5N,
		spE20W7O,
		spE20W10O,
		//spE20W10N,
		//spE25W5N,
		//spE30W5N,
		//spE40W5N,
		spE25W5O,
		spE30W5O,
		spE40W5O,
		//spAlgo,
	};
	std::cout << "Kendall Tau distance:\n";
	for (int p = 0, temSize = v_m_temp.size(); p < temSize; ++p){
		string temp_name = v_m_temp[p].first;
		double dis = .0;
		// since Kendall(t1, t1) = 0;
		for (int j = 0; j < temSize; ++j){
			dis += kendallSmallN(&v_m_temp[p].second[0],
				v_m_temp[p].first,
				&v_m_temp[j].second[0],
				v_m_temp[j].first, len, isDisplay);
		}
		// normalized tau distance.
		dis /= double(temSize - 1);
		std::cout << temp_name << " V.S. all others : "
			<< dis << std::endl;
		of << temp_name << " V.S. all others : "
			<< dis << std::endl;
	}


	of << std::endl << std::endl;
	for (int p = 0, temSize = v_m_temp.size(); p < temSize; ++p){
		string temp_name = v_m_temp[p].first;
		of << temp_name << std::endl;
	}

	of << std::endl << std::endl;
	for (int p = 0, temSize = v_m_temp.size(); p < temSize; ++p){
		string temp_name = v_m_temp[p].first;
		double dis = .0;
		// since Kendall(t1, t1) = 0;
		for (int j = 0; j < temSize; ++j){
			dis += kendallSmallN(&v_m_temp[p].second[0],
				v_m_temp[p].first,
				&v_m_temp[j].second[0],
				v_m_temp[j].first, len, false);
		}
		// normalized tau distance.
		dis /= double(temSize - 1);
		of << dis << std::endl;
	}
}
