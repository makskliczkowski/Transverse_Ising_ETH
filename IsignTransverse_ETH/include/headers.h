#pragma once

//#ifndef 
#include "config.hpp"
#include "commons.hpp"
#include "random_and_disorder/disorder.hpp"


extern int num_of_threads;													// number of threads
extern int anderson_dim;
extern std::mt19937::result_type seed_global;
extern randomGen my_gen;

using namespace std;



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
	//std::cout << "Element not found" << std::endl;
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
	NO_OVERFLOW(u64 maxPower = BinaryPowers[L - int32_t(1)];);
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

inline u64 binomial(int n, int k) {
	if (k == 0 || k == n)
		return 1;
	return binomial(n - 1, k - 1) + binomial(n - 1, k);
}
// ----------------------------------------------------------------------------- VECTORS HANDLING ----------------------------------------------------------------------------- 
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
	else if (_BC == 1) {
		iota(neis.begin(), neis.begin() + (L - int32_t(corr_len)), corr_len);
	} else
		throw "Not enough cases for me\n";
	return neis;
}

// ----------------------------------------------------------------------------- PRINTERS -----------------------------------------------------------------------------

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
		//std::cout << "couldn't open a file: " + filename << std::endl;
		return false;
	}
	return true;
}

inline void crash(bool shouldIstayorshouldIgo, std::string message = "execution terminated") {
	if (shouldIstayorshouldIgo) {
		std::cout << message;
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
		data[i] = arma::vec(num_rows, arma::fill::zeros);
	
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


template <typename ... _ty>
inline void save_to_file(std::string name, const arma::vec& x, const arma::vec& y, _ty... args) {
	if(x.size() != y.size()){
		std::cout << "Incompatible dimensions: " << x.size() << "vs.\t" << y.size() << std::endl;
		assert(false);
	}
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
	my_gen = randomGen(seed);
	std::cout << "test randoms \n" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << std::endl;
	my_gen = randomGen(seed);
	std::cout << "reset seed!" << std::endl;
	std::cout << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << std::endl;
	std::cout << "Same? Good continue!\n\n";
}
// ----------------------------------------------------------------------------- DISTRIBUTION AND DATASET RELATED FUNCTIONS -----------------------------------------------------------------------------

// ------------------------------------- definitions
arma::vec non_uniform_derivative(const arma::vec& x, const arma::vec& y);
arma::vec log_derivative(const arma::vec& x, const arma::vec& y);
double simpson_rule(double a, double b, int n, const arma::vec& f);
template<typename _type>
_type simpson_rule(const arma::vec& x, const arma::Col<_type>& f);
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
inline void setDistElem(arma::vec& dist, double min, double step, double elem) {
	const int idx = getDistIdx(min, step, elem);
	if (idx >= 0 && idx < dist.size()) dist(idx) += 1;
}

/// <summary>
/// Normalises the distribution to 1 integral
/// </summary>
inline arma::vec normalise_dist(const arma::vec& distribution, double _min, double _max) {
	return distribution / simpson_rule(_min, _max, (int)distribution.size() - 1, distribution);
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
		tmp = tmp * tmp;
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
		std::cout << "\r";															// Bring cursor to start of line
		std::cout << firstPartOfpBar;												// Print out first part of pBar
		for (int a = 0; a < amountOfFiller; a++) {								// Print out current progress
			std::cout << pBarFiller;
		}
		std::cout << pBarUpdater[currUpdateVal];
		for (int b = 0; b < pBarLength - amountOfFiller; b++) {					// Print out spaces
			std::cout << " ";
		}
		std::cout << lastPartOfpBar;												// Print out last part of progress bar
		std::cout << " (" << (int)(100 * (currentProgress / neededProgress)) << "%)";	// This just prints out the percent
		std::cout << std::flush;
		currUpdateVal += 1;
	}
	void printWithTime(const std::string& message, double percentage) {
#pragma omp critical
		{
			std::cout << "\t\t\t\t-> time: " << tim_s(timer) << message << " : \n";
			this->print();
			std::cout << std::endl;
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



#include "Lanczos/Lanczos.hpp"