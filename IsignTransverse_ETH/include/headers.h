#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <algorithm>
// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
//-------
#include <armadillo>
//#include <mkl.h>
#include <cassert> // assert terminates program
#include <omp.h>
#include <time.h>
#include <utility> // auto, etc. 
#include <memory> // smart ptr
#include <thread>
#include<queue>
#include<mutex>
#include<condition_variable>
#include<functional>
#include<future>    
#include <bitset> // binary data type
#include "random.h"

using namespace std;
using namespace arma;

typedef unsigned long long int u64;
typedef std::complex<double> cpx;
typedef std::unique_ptr<std::vector<u64>> my_uniq_ptr;

// User makros
#define im cpx(0.0,1.0)
#define out std::cout << std::setprecision(16) << std::fixed
#define num_of_threads 16

#define memory_over_performance false // optimized by size --true-- (memory usage shortage) or performance --false--

extern double pi;
extern double T; // temperature for Sq calculations
extern double dT; // temperature increment
extern double T_end; // temperature range (dT, T_end)

//----------------------------------------------------------------------------------------------
//--------------------------------------------------TOOLS---------------------------------------

/// <summary>
/// Fiunding index of base vector in mapping to reduced basis
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="arr">arary/vector conataing the mapping to the reduced basis</param> 
/// <param name="l_point">left maring for binary search</param>
/// <param name="r_point">right margin for binary search</param>
/// <param name="element">element to search in the array</param>
/// <returns></returns>
template<typename T>
inline u64 binary_search(vector<u64>& arr, u64 l_point, u64 r_point, T element) {
	if (r_point >= l_point) {
		u64 middle = l_point + (r_point - l_point) / 2;
		if (arr->at(middle) == element) return middle;
		else if (arr->at(middle) > element) return binary_search(arr, l_point, middle - 1, element);
		else return binary_search(arr, middle + 1, r_point, element);
	}
	//u64 j = findElement(arr, element);
	out << "Element " << element << " not present in the array" << endl;
	assert(false);
	return -1;
}
/// <summary>
/// Conversion to binary system
/// </summary>
/// <param name="idx">numner for conversion</param>
/// <param name="vec">vector containing the binary string</param>
inline void int_to_binary(u64 idx, std::vector<int>& vec) {
	u64 temp = idx;
	for (int k = 0; k < vec.size(); k++) {
		vec[vec.size() - 1 - k] = static_cast<int>(temp % 2);
		temp = static_cast<u64>((double)temp / 2.);
	}
}
/// <summary>
/// conversion from binary to integer
/// </summary>
/// <param name="vec">binary string</param>
/// <returns>unsigned long long integerv </returns>
inline u64 binary_to_int(vector<int>& vec) {
	u64 val = 0;
	for (int k = 0; k < vec.size(); k++) {
		val += vec[vec.size() - 1 - k] * (u64)std::pow(2, k);
	}
	return val;
}


/**
 * Overriding the ostream operator for pretty printing vectors.
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
	if (vec.size() != 0) {
		std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(os));
		os << vec.back();
	}
	return os;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, std::unique_ptr<std::vector<T>> vec) {
	if (vec->size() != 0) {
		std::copy(vec->begin(), vec->end() - 1, std::ostream_iterator<T>(os, ""));
		os << vec->back();
	}
	return os;
}
