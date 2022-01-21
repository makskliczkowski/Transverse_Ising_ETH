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
#include <unordered_map>
// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
#define ARMA_ALLOW_FAKE_GCC
#define ARMA_NO_DEBUG
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
#include <future>
#include <mutex>
//#include <condition_variable>
#include <functional>
#include <type_traits>
//#include <execution>
#ifdef __has_include
#  if __has_include(<filesystem>)
#    include <filesystem>
#    define have_filesystem 1
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
#    define have_filesystem 1
#    define experimental_filesystem
namespace fs = std::experimental::filesystem;
#  else
#    define have_filesystem 0
#  endif
#endif
#include "random.h"
#ifdef _MSC_VER
	#include <intrin.h>
	#include <nmmintrin.h>
	#define __builtin_popcount __popcnt
	#define __builtin_popcountll _mm_popcnt_u64
#endif

extern std::random_device rd;
extern std::mt19937::result_type seed;
extern std::mt19937_64 gen;
// ----------------------------------------------------------------------------- namespaces -----------------------------------------------------------------------------
using namespace std;
using namespace arma;
using clk = std::chrono::system_clock;
//namespace exec = std::execution;

// ----------------------------------------------------------------------------- definitions -----------------------------------------------------------------------------
static const char* kPathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif

using cpx = std::complex<double>;
using op_type = std::function<std::pair<cpx, u64>(u64, int, std::vector<int>)>;

template<class T>
using v_3d = std::vector<std::vector<std::vector<T>>>;											// 3d double vector
template<class T>
using v_2d = std::vector<std::vector<T>>;														// 2d double vector
template<class T>
using v_1d = std::vector<T>;																	// 1d double vector

// ----------------------------------------------------------------------------- User compiler macro -----------------------------------------------------------------------------
#if !defined(OPERATOR)
	#define OPERATOR
#endif
#if !defined(DEGENERACIES)
	//#define DEGENERACIES
#endif
#define im cpx(0.0,1.0)
#define stout std::cout << std::setprecision(8) << std::fixed								// standard outstream
#define memory_over_performance false															// optimized by size --true-- (memory usage shortage) or performance --false--

// ----------------------------------------------------------------------------- Macros to generate the lookup table (at compile-time) -----------------------------------------------------------------------------
#define R2(n) n, n + 2*64, n + 1*64, n + 3*64
#define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
#define REVERSE_BITS R6(0), R6(2), R6(1), R6(3)
#define ULLPOW(k) 1ULL << k
#define RETURNS(...) -> decltype((__VA_ARGS__)) { return (__VA_ARGS__); }
// ----------------------------------------------------------------------------- lookup table to store the reverse of each index of the table -----------------------------------------------------------------------------
// The macro `REVERSE_BITS` generates the table
const u64 lookup[256] = { REVERSE_BITS };

const v_1d<u64> BinaryPowers = { ULLPOW(0), ULLPOW(1), ULLPOW(2), ULLPOW(3),
								ULLPOW(4), ULLPOW(5), ULLPOW(6), ULLPOW(7),
								ULLPOW(8), ULLPOW(9), ULLPOW(10), ULLPOW(11),
								ULLPOW(12), ULLPOW(13), ULLPOW(14), ULLPOW(15),
								ULLPOW(16), ULLPOW(17), ULLPOW(18), ULLPOW(19),
								ULLPOW(20), ULLPOW(21), ULLPOW(22), ULLPOW(23),
								ULLPOW(24), ULLPOW(25), ULLPOW(26), ULLPOW(27),
								ULLPOW(28), ULLPOW(29), ULLPOW(30), ULLPOW(31) }; // vector containing powers of 2 from 2^0 to 2^(L-1)

// ----------------------------------------------------------------------------- CONSTANTS -----------------------------------------------------------------------------
extern int num_of_threads;													// number of threads
constexpr long double pi = 3.141592653589793238462643383279502884L;			// it is me, pi
constexpr long double two_pi = 2 * 3.141592653589793238462643383279502884L;	// it is me, 2pi
const auto global_seed = std::random_device{}();							// global seed for classes
const std::string kPSep = std::string(kPathSeparator);
// ----------------------------------------------------------------------------- TIME FUNCTIONS -----------------------------------------------------------------------------

inline double tim_s(clk::time_point start) {
	return double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(\
		std::chrono::system_clock::now() - start)).count()) / 1000.0;
}
// ----------------------------------------------------------------------------- TOOLS -----------------------------------------------------------------------------
enum class coordinate {
	x, 
	y,
	z
};

/// <summary>
/// Calculates the sign of a value
/// </summary>
template <typename T> int sgn(T val) {
	return int(T(0) < val) - int(val < T(0));
}
// ----------------------------------------------------------------------------- STRING BASED TOOLS DECLARATIONS -----------------------------------------------------------------------------
// weird: interfase to enforce same types in variadic templates
template <typename _baseTy>
struct types {
	template<typename... _otherTy>
	using convertible = typename std::enable_if<std::conjunction<std::is_convertible<_otherTy, _baseTy>...>::value>::type;
};

// ---------------------------------- definitions
bool isNumber(const string& str);

std::vector<std::string> split_str(std::string s, std::string delimiter);

// ---------------------------------- templates

/// <summary>
/// Changes a value to a string with a given precison
/// </summary>
/// <param name="n">number of decimal places</param>
/// <returns>string of a number with given precision</returns>
template <typename T>
std::string to_string_prec(const T a_value, const int n = 3) {
	std::ostringstream outie;
	outie.precision(n);
	outie << std::fixed << a_value;
	return outie.str();
}

// ----------------------------------------------------------------------------- BINARY TOOLS -----------------------------------------------------------------------------

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

/// Template instance

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

// ---------------------------------- inlines

/// <summary>
/// Rotates the binary representation of the input decimal number by one left shift
/// </summary>
/// <param name="n"> number to rotate </param>
/// <param name="maxPower"> maximal power of 2 </param>
/// <returns> rotated number </returns>
inline u64 rotate_left(u64 n, int L) {
	u64 maxPower = BinaryPowers[L - 1];
	return (n >= maxPower) ? (((int64_t)n - (int64_t)maxPower) * 2 + 1) : n * 2;
}

/// <summary>
/// Check the k'th bit
/// </summary>
/// <param name="n">Number on which the bit shall be checked</param>
/// <param name="k">number of bit (from 0 to 63)</param>
/// <returns>Bool on if the bit is set or not</returns>
inline bool checkBit(u64 n, int k) {
	return n & (1ULL << k);
}

/// <summary>
/// flip the bits in the number. The flipping is done via substracting the maximal number we can get for a given bitnumber
/// </summary>
/// <param name="n">number to be flipped</param>
/// <param name="maxBinaryNum">maximal power of 2 for given bit number(maximal length is 64 for ULL)</param>
/// <returns>flipped number</returns>
inline u64 flip(u64 n, int L) {
	return BinaryPowers[L] - n - 1;
}

/// <summary>
/// Flip the bit on k'th site and return the number it belongs to. The bit is checked from right to left!
/// </summary>
/// <param name="n">number to be checked</param>
/// <param name="kthPower">precalculated power of 2 for k'th site</param>
/// <param name="k">k'th site for flip to be checked</param>
/// <returns>number with k'th bit from the right flipped</returns>
inline u64 flip(u64 n, u64 kthPower, int k) {
	return checkBit(n, k) ? (int64_t(n) - (int64_t)kthPower) : (n + kthPower);
}

/// <summary>
/// Function that calculates the bit reverse, note that 64 bit representation
/// is now taken and one has to be sure that it doesn't exceede it (which it doesn't, we sure)
/// </summary>
/// <param name="L">We need to know how many bits does the number really take because the function can take up to 64</param>
/// <returns>number with reversed bits moved to be maximally of size L again</returns>
inline u64 reverseBits(u64 n, int L) {
	u64 rev = (lookup[n & 0xffULL] << 56) |					// consider the first 8 bits
		(lookup[(n >> 8) & 0xffULL] << 48) |				// consider the next 8 bits
		(lookup[(n >> 16) & 0xffULL] << 40) |				// consider the next 8 bits
		(lookup[(n >> 24) & 0xffULL] << 32) |				// consider the next 8 bits
		(lookup[(n >> 32) & 0xffULL] << 24) |				// consider the next 8 bits
		(lookup[(n >> 40) & 0xffULL] << 16) |				// consider the next 8 bits
		(lookup[(n >> 48) & 0xffULL] << 8) |				// consider the next 8 bits
		(lookup[(n >> 54) & 0xffULL]);						// consider last 8 bits
	return (rev >> (64 - L));								// get back to the original maximal number
}

inline std::function<u64(u64, int)> multiply_operators(const std::function<u64(u64, int)>& A, const std::function<u64(u64, int)>& B) {
	std::function<u64(u64, int)> result = [A, B](u64 n, int L) { return A(B(n, L), L); };
	return result;
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
		vec[size - 1 - k] = temp % 2;
		temp = temp / 2.;
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
/// change vector of bools to int using the vector of powers precalculated
/// </summary>
inline u64 binary_to_int(const vector<bool>& vec, const v_1d<u64>& powers) {
	u64 val = 0;
	const u64 size = vec.size();
	for (int k = 0; k < size; k++)
		val += static_cast<u64>(vec[size - 1 - k]) * powers[k];
	return val;
}

// ----------------------------------------------------------------------------- VECTORS HANDLING -----------------------------------------------------------------------------

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

/// <summary>
/// Creates a random vector of custom length using the random library and the merson-twister (?) engine
/// </summary>
/// <param name="N"> length of the generated random vector </param>
/// <returns> returns the custom-length random vector </returns>
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
inline std::vector<int> get_neigh_vector(int _BC, int L, int corr_len) {
	v_1d<int> neis(L, -1);
	if (_BC == 0) {
		iota(neis.begin(), neis.end(), 0);
		std::rotate(neis.begin(), neis.begin() + corr_len, neis.end());
	}
	else if (_BC == 1)
		iota(neis.begin(), neis.begin() + (L - corr_len), corr_len);
	else
		throw "Not enough cases for me\n";
	return neis;
}

// ----------------------------------------------------------------------------- PRINTERS -----------------------------------------------------------------------------

/// <summary>
/// Overriding the ostream operator for pretty printing vectors.
/// </summary>
/// <typeparam name="T"> writing out </typeparam>
/// <param name="os"> designed outstream </param>
/// <param name="vec"> vector variable to print </param>
/// <returns></returns>
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
	int counter = 0;
	if (vec.size() != 0) {
		std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(os, "\t\t"));
		os << vec.back() << ' ';
		//for (int i = 0; i < vec.size(); i++) {
		//	os << vec[i] << "\t\t";
		//	//counter++;
		//	//if (counter % 8 == 0) {
		//	//	os << "\t\t";
		//	//	counter = 0;
		//	//}
		//}
	}
	else
		os << "Empty container!" << endl;
	return os;
}

/// <summary>
///
/// </summary>
/// <param name="filename"></param>
/// <param name="mode"></param>
/// <returns></returns>
template <typename T>
inline bool openFile(T& file, std::string filename, std::ios_base::openmode mode = std::ios::out) {
	file.open(filename, mode);
	if (!file.is_open()) {
		stout << "couldn't open a file: " + filename + "\n";
		return false;
	}
	return true;
}

inline void crash(bool shouldIstayorshouldIgo, std::string message = "execution terminated") {
	if (shouldIstayorshouldIgo) {
		stout << message;
		exit(1);
	}
}
/// <summary>
///
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="output"></param>
/// <param name="elements"></param>
/// <param name="seperator"></param>
template <typename _ty>
inline void print(std::ostream& output, std::string separator, arma::u16 width, _ty arg) {
	output.width(width); output << arg << separator;
}
template <typename _T, typename... _ty>
inline void print(std::ostream& output, std::string separator, arma::u16 width, _T arg, _ty... rest) {
	print(output, separator, width, arg);
	print(output, separator, width, rest...);
}
template <typename... _ty>
inline void printSeparated(std::ostream& output, std::string separator, arma::u16 width, bool endline, _ty... args) {
	print(output, separator, width, args...);
	if (endline) output << std::endl;
}

inline auto readFromFile(std::ifstream& input, std::string filename) {
	std::string datarow;
	std::vector<arma::vec> data;
	// find dimensionality of data
	int num_cols = 0;
	int num_rows = 0;
	bool isopen = openFile(input, filename, ios::in);
	if (!isopen) return std::vector<arma::vec>();
	while (std::getline(input, datarow)) {
		std::istringstream ss(datarow);
		if (num_rows == 0) {
			double value;
			while (ss >> value)	num_cols++;
			data = std::vector<arma::vec>(num_cols);
		}
		num_rows++;
	}
	for (int i = 0; i < num_cols; i++)
		data[i] = arma::vec(num_rows, fill::zeros);
	
	// load data
	input.close();
	openFile(input, filename, ios::in);
	int j = 0;
	while(std::getline(input, datarow)){
		std::istringstream ss(datarow);
		int i = -1;
		double value;
		while (ss >> value) {
			data[++i](j) = value;
		}
		j++;
	}
	return data;
}

inline void createDirs(const std::string& dir) {
	if (!fs::is_directory(dir) || !fs::exists(dir))
		fs::create_directories(dir);
}
template <typename... _Ty>
inline void createDirs(const std::string& dir, const _Ty&... dirs) {
	createDirs(dir);
	createDirs(dirs...);
}


template <typename ... _ty>
inline void save_to_file(std::string name, const arma::vec& x, const arma::vec& y, _ty... args) {
	assert(x.size() == y.size() && "Incompatible dimensions");
	std::ofstream file;
	openFile(file, name, ios::out);
	for (int i = 0; i < x.size(); i++) {
		if (i == 0) printSeparated(file, "\t", 12, true, x(i), y(i), args...);
		else		printSeparated(file, "\t", 12, true, x(i), y(i));
	}
	file.close();
}
template <typename ... _ty>
inline void save_to_file(std::string name, const arma::vec& x, const arma::vec& y, const arma::vec& z, _ty... args) {
	assert(((x.size() == y.size()) && (x.size() == z.size())) && "Incompatible dimensions");
	std::ofstream file;
	openFile(file, name, ios::out);
	for (int i = 0; i < x.size(); i++) {
		if (i == 0) printSeparated(file, "\t", 12, true, x(i), y(i), z(i), args...);
		else		printSeparated(file, "\t", 12, true, x(i), y(i), z(i));
	}
	file.close();
}
// ----------------------------------------------------------------------------- MAKS' IDEAS -----------------------------------------------------------------------------
/// <summary>
/// Sorts the vector and saves the permutation with a lambda like function compare
/// </summary>
template <typename T, typename Compare> inline std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) {
			return compare(vec[i], vec[j]);
		});
	return p;
}

/// <summary>
/// Applies permutation on a given vector
/// </summary>
template <typename T> inline void apply_permutation(std::vector<T>& vec, const std::vector<std::size_t>& p) {
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

template <typename T>
inline T variance(T value, T average, int norm) {
	return std::sqrt((value / norm - average * average) / norm);
}

template <typename _Ty>
inline _Ty matrixVariance(const arma::Mat<_Ty>& mat) {
	_Ty var = 0, mean = 0;
#pragma omp parallel for reduction(+: var, mean) collapse(2)
	for (long int n = 0; n < mat.n_cols; n++)
		for (long int m = 0; m < mat.n_rows; m++) {
			var += mat(n, m) * mat(n, m);
			mean += mat(n, m);
		}
	var /= double(mat.n_cols*mat.n_rows);
	mean /= double(mat.n_cols*mat.n_rows);
	return var - mean * mean;
}
inline void checkRandom(unsigned int seed) {
	gen = std::mt19937_64(seed);
	std::uniform_int_distribution<int> dist;
	stout << "test randoms \n" << dist(gen) << "\t" << dist(gen) << "\t" << dist(gen) << std::endl;
	gen = std::mt19937_64(seed);
	stout << "reset seed!" << std::endl;
	stout << dist(gen) << "\t" << dist(gen) << "\t" << dist(gen) << std::endl;
	stout << "Same? Good continue!\n\n";
}
// ----------------------------------------------------------------------------- DISTRIBUTION AND DATASET RELATED FUNCTIONS -----------------------------------------------------------------------------

// ------------------------------------- definitions
double simpson_rule(double a, double b, int n, const arma::vec& f);
double simpson_rule(const arma::vec& x, const arma::vec& f);
double binder_cumulant(const arma::vec& arr_in);														// calculate binder cumulant of dataset
arma::vec get_NonDegenerated_Elements(const arma::vec& arr_in);											// compute non-unique values in dataset
// ------------------------------------- inlines

/// <summary>
/// From minimum and a step finds an index in a given probability distribution
/// </summary>
inline int getDistIdx(double min, double step, double elem) {
	return static_cast<int>((elem + abs(min)) / step);
}

/// <summary>
/// Sets up the element in a given distribution dist
/// </summary>
inline void setDistElem(vec& dist, double min, double step, double elem) {
	const int idx = getDistIdx(min, step, elem);
	if (idx >= 0 && idx < dist.size()) dist(idx) += 1;
}

/// <summary>
/// Normalises the distribution to 1 integral
/// </summary>
inline arma::vec normalise_dist(const arma::vec& distribution, double _min, double _max) {
	return distribution / simpson_rule(_min, _max, distribution.size() - 1, distribution);
}

/// <summary>
/// Calculates the i'th moment
/// </summary>
/// <param name="power">number of the moment</param>
inline double moment(const arma::vec& arr_in, int power) {
	return arma::mean(arma::pow(arr_in - mean(arr_in), power));
}

/// <summary>
/// Calculates the kurtosis of a given vector
/// </summary>
/// <param name="arr_in"></param>
/// <returns></returns>
inline double kurtosis(const arma::vec& arr_in) {
	double dev_inv = 1.0 / std::pow(arma::stddev(arr_in), 4);
	return moment(arr_in, 4) * dev_inv - 3.0;
}

/// <summary>
/// Kurtosis in different way using openMP
/// </summary>
inline double kurtosis_diff(const arma::vec& arr_in) {
	double mean = arma::mean(arr_in);
	double fourth = 0;
	double second = 0;
	double counter = 0;
#pragma omp parallel for reduction (+: fourth, second, counter)
	for (int i = 0; i < arr_in.size(); i++) {
		double tmp = arr_in(i) - mean;
		tmp *= tmp;
		fourth += tmp * tmp;
		second += tmp;
		counter += 1;
	}
	fourth = fourth / counter;
	second = second / counter;
	return fourth / (second * second) - 3.0;
}

/// <summary>
/// Compute Gaussina function at value x with input mean and variance
/// </summary>
/// <param name="x"> input value </param>
/// <param name="meam"> mean value of gaussian function </param>
/// <param name="std_dev"> standard deviation of gaussian </param>
template <class T>
inline T gaussian(T x, double mean, double std_dev) {
	T exponent = (x - mean) / std_dev;
	return 1.0 / (std::sqrt(two_pi) * std_dev) * exp(-pow(exponent, 2) / 2.0);
}

// --------------------------------- PROGRESS BAR ---------------------------------
class pBar {
public:
	void update(double newProgress) {
		currentProgress += newProgress;
		amountOfFiller = (int)((currentProgress / neededProgress) * (double)pBarLength);
	}
	void print() {
		currUpdateVal %= pBarUpdater.length();
		stout << "\r";															// Bring cursor to start of line
		stout << firstPartOfpBar;												// Print out first part of pBar
		for (int a = 0; a < amountOfFiller; a++) {								// Print out current progress
			stout << pBarFiller;
		}
		stout << pBarUpdater[currUpdateVal];
		for (int b = 0; b < pBarLength - amountOfFiller; b++) {					// Print out spaces
			stout << " ";
		}
		stout << lastPartOfpBar;												// Print out last part of progress bar
		stout << " (" << (int)(100 * (currentProgress / neededProgress)) << "%)";	// This just prints out the percent
		stout << std::flush;
		currUpdateVal += 1;
	}
	void printWithTime(const std::string& message, double percentage) {
#pragma omp critical
		{
			stout << "\t\t\t\t-> time: " << tim_s(timer) << message << " : \n";
			this->print();
			stout << std::endl;
		}
		this->update(percentage);
	}
	// constructor
	pBar() {
		timer = std::chrono::system_clock::now();
		amountOfFiller = 0;
	}
private:
	// --------------------------- STRING ENDS
	std::string firstPartOfpBar = "\t\t\t\t[";
	std::string lastPartOfpBar = "]";
	std::string pBarFiller = "|";
	std::string pBarUpdater = "/-\\|";
	// --------------------------- PROGRESS
	clk::time_point timer;														// inner clock
	int amountOfFiller;															// length of filled elements
	int pBarLength = 50;														// length of a progress bar
	int currUpdateVal = 0;														//
	double currentProgress = 0;													// current progress
	double neededProgress = 100;												// final progress
};