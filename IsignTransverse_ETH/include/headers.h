#pragma once
#include "config.hpp"
#include "commons.h"
#include "random.h"
#include "digamma.h"

extern int num_of_threads;													// number of threads
extern int anderson_dim;
extern std::mt19937::result_type seed_global;
extern randomGen my_gen;

// ----------------------------------------------------------------------------- namespaces -----------------------------------------------------------------------------
using namespace std;
//namespace exec = std::execution;



// ----------------------------------------------------------------------------- User compiler macro -----------------------------------------------------------------------------


/// <summary>
/// Calculates the sign of a value
/// </summary>
template <typename T> int sgn(T val) {
	return int(T(0) < val) - int(val < T(0));
}
// ----------------------------------------------------------------------------- STRING BASED TOOLS DECLARATIONS -----------------------------------------------------------------------------
// weird: interfase to enforce same types in variadic templates

// ---------------------------------- definitions
bool isNumber(const string& str);

std::vector<std::string> split_str(std::string s, std::string delimiter);

// ---------------------------------- templates

//<! finds the order of magnitude of number +1 (only for <1 numbers to find filename format)
template <typename T>
inline
int order_of_magnitude(const T a_value) {
	
	if(a_value < 1.0 && a_value != 0){
		T m = std::abs(std::log10(std::abs(a_value)));
		return int(std::max(std::ceil(m) + 1., 2.));
	}
	else return 2;
}
/// <summary>
/// Changes a value to a string with a given precison
/// </summary>
/// <param name="n">number of decimal places</param>
/// <returns>string of a number with given precision</returns>
template <typename T>
inline
std::string to_string_prec(const T a_value, int n = -1) {
	if(n < 0)
		n = order_of_magnitude(a_value);
	std::cout << n << std::endl;
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
		//stout << "couldn't open a file: " + filename << std::endl;
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
// ----------------------------------------------------------------------------- MAKS' IDEAS -----------------------------------------------------------------------------
/// <summary>
/// Sorts the vector and saves the permutation with a lambda like function compare
/// </summary>
template <typename T, typename Compare> 
inline std::vector<std::size_t> sort_permutation(const T& vec, Compare compare) {
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
template <typename T> 
inline void apply_permutation(T& vec, const std::vector<std::size_t>& p) {
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
	my_gen = randomGen(seed);
	stout << "test randoms \n" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << std::endl;
	my_gen = randomGen(seed);
	stout << "reset seed!" << std::endl;
	stout << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << "\t" << my_gen.random_uni<double>(0., 1.) << std::endl;
	stout << "Same? Good continue!\n\n";
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


//! ----------------------------------------------------------------------------- ARMADILLO HELPERS -----------------------------------------------------------------------------
//<! calculate commutator of two input matrix types, which have overloaded * operator
inline std::string matrix_size(u64 dim){
	 if(dim < 1e3)
	 	return std::to_string(dim) + " bytes";
	 else if(dim < 1e6)
	 	return to_string_prec(dim / 1e3, 2) + " kB";
	 else if(dim < 1e9)
	 	return to_string_prec(dim / 1e6, 2) + " MB";
	 else if(dim < 1e12)
	 	return to_string_prec(dim / 1e9, 2) + " GB";
	else 
	 	return to_string_prec(dim / 1e12, 2) + " TB";
}


template <typename matrix>
matrix commutator(const matrix& A, const matrix& B)
	{ return A * B - B * A; }

//! -------------------------------------------------------- cast non-cpx to cpx types
template <typename _ty>
arma::Col<std::complex<_ty>> cpx_real_vec(const arma::Col<_ty>& input){ 
	size_t size = input.size();
	return arma::Col<std::complex<_ty>>(input, arma::Col<_ty>(size, arma::fill::zeros));
}
template <typename _ty>
arma::Col<std::complex<_ty>> cpx_imag_vec(const arma::Col<_ty>& input) {
	size_t size = input.size();
	return arma::Col<std::complex<_ty>>(arma::Col<_ty>(size, arma::fill::zeros), input);
}
template <typename _ty>
arma::Col<std::complex<_ty>> cpx_real_vec(const arma::subview_col<_ty>& input) {
	size_t size = input.n_elem;
	return arma::Col<std::complex<_ty>>(input, arma::Col<_ty>(size, arma::fill::zeros));
}
template <typename _ty>
arma::Col<std::complex<_ty>> cpx_imag_vec(const arma::subview_col<_ty>& input) {
	size_t size = input.n_elem;
	return arma::Col<std::complex<_ty>>(arma::Col<_ty>(size, arma::fill::zeros), input);
}


template <typename _type>
inline
arma::cx_vec cast_cx_vec(const arma::Col<_type>& state);

template <>
inline arma::cx_vec cast_cx_vec(const arma::vec& state)
	{ return cpx_real_vec(state); }
template <>
inline arma::cx_vec cast_cx_vec(const arma::cx_vec& state)
	{ return state; }
//! -------------------------------------------------------- dot product for different input types (cpx and non-cpx)

 template <typename _ty, 
	 template <typename> class _COLVEC1,
	 template <typename> class _COLVEC2 
 >
_ty dot_prod(const _COLVEC1<_ty>& left, const _COLVEC2<_ty>& right)
{ 
	 static_assert(_traits::is_any_of_v<_COLVEC1<_ty>, arma::Col<_ty>, arma::subview_col<_ty>>
		 && _traits::is_any_of_v<_COLVEC2<_ty>, arma::Col<_ty>, arma::subview_col<_ty>>,
		 "Dot product only valid for arma::Col and arma::subview classes");
	return arma::cdot(left, right); 
}																			
																											
template <typename _ty,
	template <typename> class _COLVEC1,
	template <typename> class _COLVEC2
>
std::complex<_ty> dot_prod(const _COLVEC1<_ty>& left, const _COLVEC2<std::complex<_ty>>& right)
{
	static_assert(_traits::is_any_of_v<_COLVEC1<_ty>, arma::Col<_ty>, arma::subview_col<_ty>>
		&& _traits::is_any_of_v<_COLVEC2<std::complex<_ty>>, arma::Col<std::complex<_ty>>, arma::subview_col<std::complex<_ty>>>,
		"Dot product only valid for arma::Col and arma::subview classes");
	return arma::cdot(cpx_real_vec(left), right);
}
																											
template <typename _ty,
	template <typename> class _COLVEC1,
	template <typename> class _COLVEC2
>
std::complex<_ty> dot_prod(const _COLVEC1<std::complex<_ty>> & left, const _COLVEC2<_ty> & right)
{
	static_assert(_traits::is_any_of_v<_COLVEC1<std::complex<_ty>>, arma::Col<std::complex<_ty>>, arma::subview_col<std::complex<_ty>>>
		&& _traits::is_any_of_v<_COLVEC2<_ty>, arma::Col<_ty>, arma::subview_col<_ty>>,
		"Dot product only valid for arma::Col and arma::subview classes"); 
	return arma::cdot(left, cpx_real_vec(right));
}															
											
template <typename _ty>
inline arma::Col<_ty> exctract_vector_between_values(
	const arma::Col<_ty>& input_vec,	//<! input vector to exctract data from (assumed sorted)
	_ty start, 							//<! first value of new vector (if lower than lowest in input_vec than taking from beggining)
	_ty end								//<! last element to copy data
) {
	arma::Col<_ty> output;
	for (auto& it : input_vec) {
		if (it >= start && it <= end) {
			int size = output.size();
			output.resize(size + 1);
			output(size) = it;
		}
	}
	return output;
}
inline arma::vec exctract_vector(
	const arma::vec& input_vec,	//<! input vector to exctract data from (assumed sorted)
	u64 start, 	//<! first value of new vector (if lower than lowest in input_vec than taking from beggining)
	u64 end		//<! last element to copy data
) {
	arma::vec output(end - start);
#pragma omp parallel for
	for (int k = start; k < end; k++) 
		output(k - start) = input_vec(k);
	return output;
}

template <typename _type>
inline
arma::sp_cx_mat cast_cx_sparse(const arma::SpMat<_type>& mat);

template <>
inline arma::sp_cx_mat cast_cx_sparse(const arma::sp_mat& mat)
{
	arma::sp_cx_mat ret(mat.n_rows, mat.n_cols);
	ret.set_real(mat);
	return ret;
}
template <>
inline arma::sp_cx_mat cast_cx_sparse(const arma::sp_cx_mat& mat)
	{ return mat; }
//general_dot_prod(arma::Col,			arma::Col		 );
//general_dot_prod(arma::subview_col, arma::Col		 );
//general_dot_prod(arma::Col,			arma::subview_col);
//general_dot_prod(arma::subview_col, arma::subview_col);

#include "Lanczos/Lanczos.hpp"