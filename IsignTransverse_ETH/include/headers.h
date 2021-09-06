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
#include <filesystem>
// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
#define ARMA_ALLOW_FAKE_GCC
//-------
#include <armadillo>
#include <iterator>
//#include <mkl.h>
#include <cassert> // assert terminates program
#include <omp.h>
#include <ctime>
#include <utility> // auto, etc.
#include <memory> // smart ptr
#include <thread>
#include <mutex>
//#include <condition_variable>
#include <functional>
#include <execution>

#include "random.h"




// - - - - - namespaces - - - - -
using namespace std;
using namespace arma;
namespace fs = std::filesystem;
namespace exec = std::execution;

//- - - - - definitions - - - - -
static const char* kPathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif

using u64 = unsigned long long;
using cpx = std::complex<double>;

template<class T>
using v_3d = std::vector<std::vector<std::vector<T>>>;						// 3d double vector
template<class T>
using v_2d = std::vector<std::vector<T>>;									// 2d double vector
template<class T>
using v_1d = std::vector<T>;												// 1d double vector

// User compiler macro
#define im cpx(0.0,1.0)
#define stout std::cout << std::setprecision(16) << std::fixed				// standard outstream
#define memory_over_performance false										// optimized by size --true-- (memory usage shortage) or performance --false--

// Constants
extern int num_of_threads;													// number of threads
constexpr long double pi = 3.141592653589793238462643383279502884L;			// it is me, pi
constexpr long double two_pi = 2 * 3.141592653589793238462643383279502884L;	// it is me, 2pi
const auto global_seed = std::random_device{}();							// global seed for classes 

//static random_num* rn; // random number class instance
extern std::random_device rd;
extern std::mt19937::result_type seed;
extern std::mt19937_64 gen;

//--------------------------------------------------TOOLS--------------------------------------------------
template <typename T> int sgn(T val) {
	return int(T(0) < val) - int(val < T(0));
}

/* STRING BASED TOOLS DECLARATIONS */
bool isNumber(const string& str);

std::vector<std::string> split_str(std::string s, std::string delimiter);

template <typename T>
std::string to_string_prec(const T a_value, const int n = 3) {
	std::ostringstream outie;
	outie.precision(n);
	outie << std::fixed << a_value;
	return outie.str();
}

/* DEFINITIONS */
/// <summary>
/// Fiunding index of base vector in mapping to reduced basis
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="arr"> arary/vector conataing the mapping to the reduced basis </param>
/// <param name="l_point"> left maring for binary search </param>
/// <param name="r_point"> right margin for binary search </param>
/// <param name="element"> element to search in the array </param>
/// <returns></returns>
template <class T>
inline u64 binary_search(const std::vector<T>& arr, u64 l_point, u64 r_point, T element) {
	if (l_point < 0) assert(false && "What?");
	if (r_point >= arr.size()) {
		return -1;
	}
	if (r_point >= l_point) {
		u64 middle = l_point + (r_point - l_point) / 2;
		if (arr[middle] == element) return middle;
		else if (arr[middle] < element) return binary_search(arr, middle + 1, r_point, element);
		else return binary_search(arr, l_point, middle - 1, element);
	}
	return -1;
}
template <>
inline u64 binary_search(const std::vector<double>& arr, u64 l_point, u64 r_point, double element) {
	if (l_point < 0) assert(false && "What?");
	if (r_point >= arr.size()) {
		return -1;
	}
	if (r_point >= l_point) {
		u64 middle = l_point + (r_point - l_point) / 2;
		if (abs(arr[middle] - element) < 1e-12) return middle;
		else if (arr[middle] < element) return binary_search(arr, middle + 1, r_point, element);
		else return binary_search(arr, l_point, middle - 1, element);
	}
	return -1;
}

/// <summary>
/// Conversion to binary system
/// </summary>
/// <param name="idx"> numner for conversion </param>
/// <param name="vec"> vector containing the binary string </param>
inline void int_to_binary(u64 idx, std::vector<bool>& vec) {
	u64 temp = idx;
	const u64 size = vec.size();
	for (int k = 0; k < size; k++) {
		vec[size - 1 - k] = static_cast<bool>(temp % 2);
		temp = static_cast<u64>((double)temp / 2.);
	}
}

/// <summary>
/// conversion from binary to integer
/// </summary>
/// <param name="vec"> binary string </param>
/// <returns> unsigned long long integer </returns>
inline u64 binary_to_int(const vector<bool>& vec) {
	u64 val = 0;
	u64 exp = 1;
	const u64 size = vec.size();
	for (int k = 0; k < size; k++) {
		val += static_cast<u64>(vec[size - 1 - k]) * exp;
		exp *= 2;
	}
	return val;
}

/// <summary>
/// Creates a random vector of custom length using the random library and the merson-twister (?) engine
/// </summary>
/// <param name="N"> length of the generated random vector </param>
/// <returns> returns the custom-length random vector </returns>
inline vec create_random_vec(u64 N, double h = 1.0) {
	vec random_vec(N, fill::zeros);
	std::uniform_real_distribution<double> distribute(-h, h);
	// create random vector from middle to always append new disorder at lattice endpoint
	for (u64 j = 0; j <= N / 2.; j++) {
		u64 idx = N / 2. - j;
		random_vec(idx) = distribute(gen);
		idx += 2 * j;
		if (idx < N) random_vec(idx) = distribute(gen);
	}
	return random_vec;
}
inline std::vector<double> create_random_vec_std(u64 N) {
	std::vector<double> random_vec(N, 0);
	std::uniform_real_distribution<double> distribute(-1.0, 1.0);
	for (u64 j = 0; j < N; j++) {
		random_vec[j] = distribute(gen);
	}
	return random_vec;
}

/// <summary>
/// Calculate the vector that consists of a given site corr_len away for a lattice site provided by the user
/// </summary>
/// <param name="_BC">boundary conditions, 0 - PBC, 1 - OBC</param>
/// <param name="L">chain length</param>
/// <param name="corr_len">correlation length</param>
/// <returns>vector of correlation places</returns>
inline std::vector<int> get_neigh_vector(int _BC, int L, int corr_len){
	v_1d<int> neis(L, -1);
	if(_BC == 0){
		iota(neis.begin(), neis.end(), 0);
		std::rotate(neis.begin(), neis.begin() + corr_len,neis.end());
	}
	else if(_BC == 1)
		iota(neis.begin(), neis.begin() + (L-corr_len), corr_len);
	else
		throw "Not enough cases for me\n";

}

/// <summary>
/// find non-unique elements in input array and store only elemetns, which did not have duplicates
/// </summary>
/// <param name="arr_in"> degenerated input array </param>
/// <returns> indices to unique elements </returns>
inline arma::vec get_NonDegenerated_Elements(const arma::vec& arr_in) {
	std::vector<double> arr_degen;
	u64 N = arr_in.size();
	for (int k = 0; k < N - 1; k++)
		if (abs(arr_in(k + 1) - arr_in(k)) < 1e-12)
			arr_degen.push_back(arr_in(k));
	if (abs(arr_in(N - 1) - arr_in(N - 2)) < 1e-12)
		arr_degen.push_back(arr_in(N - 1));

	std::vector<double> arr_unique = arma::conv_to<std::vector<double>>::from(arr_in);
	for (auto& E : arr_degen){
		auto new_end = remove_if(arr_unique.begin(), arr_unique.end(), [&](double a) {
		return abs(a - E) < 1e-12;
			});
		arr_unique.erase(new_end, arr_unique.end());
	}
	return arma::unique((vec)arr_unique);
}

/// <summary>
/// Overriding the ostream operator for pretty printing vectors.
/// </summary>
/// <typeparam name="T"> writing out </typeparam>
/// <param name="os"> designed outstream </param>
/// <param name="vec"> vector variable to print </param>
/// <returns></returns>
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
	if (vec.size() != 0) {
		std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(os, " "));
		os << vec.back() << ' ';
	}
	else
		os << "Empty container!" << endl;
	return os;
}

template <typename T>
inline std::vector<T> operator+(const std::vector<T>& A, const std::vector<T>& B) {
	std::vector<T> res(A.size());
	assert(A.size() == B.size() && "stupid cunt, wonrg size of vectors\n");
	for (int k = 0; k < A.size(); k++)
		res[k] = A[k] + B[k];
	return res;
}

template <typename T>
inline std::vector<T> operator-(const std::vector<T>& A, const std::vector<T>& B) {
	std::vector<T> res(A.size());
	assert(A.size() == B.size() && "stupid cunt, wonrg size of vectors\n");
	for (int k = 0; k < A.size(); k++)
		res[k] = A[k] - B[k];
	return res;
}

template <typename T, typename T2>
inline std::vector<T> operator*(const T2 B, const std::vector<T>& A) {
	std::vector<T> res(A.size());
	for (int k = 0; k < A.size(); k++)
		res[k] = B * A[k];
	return res;
}

/* MAKS' USELESS IDEAS */
template <typename T, typename Compare>
inline std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) {
			return compare(vec[i], vec[j]);
		});
	return p;
}

template <typename T>
inline void apply_permutation(std::vector<T>& vec, const std::vector<std::size_t>& p) {
	std::vector<bool> done(vec.size());
	for (std::size_t i = 0; i < vec.size(); ++i) {
		if (done[i]) continue;
		done[i] = true;
		std::size_t prev_j = i;
		std::size_t j = p[i];
		while (i != j) {
			std::swap(vec[prev_j], vec[j]);
			done[j] = true;
			prev_j = j;
			j = p[j];
		}
	}
}