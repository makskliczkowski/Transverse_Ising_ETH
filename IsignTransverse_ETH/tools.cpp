#include "include/headers.h"
int num_of_threads = 16;

std::random_device rd;
std::mt19937::result_type seed = static_cast<long unsigned int>(time(0)); // 87178291199L; // set constant to maintain same disorder for different sizes etc
// rd() ^ (
//	(std::mt19937::result_type)
//	std::chrono::duration_cast<std::chrono::seconds>(
//		std::chrono::system_clock::now().time_since_epoch()
//		).count() +
//	(std::mt19937::result_type)
//	std::chrono::duration_cast<std::chrono::microseconds>(
//		std::chrono::high_resolution_clock::now().time_since_epoch()
//		).count());
std::mt19937_64 gen(seed);

/* STRING BASED TOOLS */
/// <summary>
/// Checking if given string is a number
/// </summary>
/// <param name="str"> string to be checked </param>
/// <returns> true if is a number </returns>
bool isNumber(const string& str)
{
	bool found_dot = false;
	bool found_minus = false;
	for (char const& c : str) {
		if (c == '.' && !found_dot) {
			found_dot = true;
			continue;
		}
		if (c == '-' && !found_minus) {
			found_minus = true;
			continue;
		}
		if (std::isdigit(c) == 0) return false;
	}
	return true;
}
/// <summary>
/// Splits string by a given delimiter
/// </summary>
/// <param name="s"> string to be split </param>
/// <param name="delimiter"> a splitting delimiter </param>
/// <returns></returns>
std::vector<std::string> split_str(std::string s, std::string delimiter)
{
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	string token;
	vector<string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(s.substr(pos_start));
	return res;
}



// PROBABILITY BASED TOOLS
double simpson_rule(double a, double b, int n, const arma::vec& f) {
	double h = (b - a) / n;

	// Internal sample points, there should be n - 1 of them
	double sum_odds = 0.0;
#pragma omp parallel for reduction(+: sum_odds)
	for (int i = 1; i < n; i += 2) {
		int idx = int(((a + i * h) + abs(a)) / h);
		sum_odds += (idx < f.size()) ? f(idx) : 0.0;
	}

	double sum_evens = 0.0;
#pragma omp parallel for reduction(+: sum_evens)
	for (int i = 2; i < n; i += 2) {
		int idx = int(((a + i * h) + abs(a)) / h);
		sum_evens += (idx < f.size()) ? f(idx) : 0.0;
	}

	return (f(0) + f(f.size() - 1) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

DISABLE_WARNING_PUSH
DISABLE_OVERFLOW // VS19 error with simple implicit cast, instead produces warning: might overflow
arma::vec non_uniform_derivative(const arma::vec& x, const arma::vec& y) {
	const int size = y.size();
	arma::vec output(size - 1, arma::fill::zeros);
	for (int j = 1; j < size; j++) 
		output(j - 1) = (y(j) - y(j - 1)) / (x(j) - x(j - 1));
	return output;
}
arma::vec log_derivative(const arma::vec& x, const arma::vec& y){
	const int size = y.size();
	arma::vec output(size - 1, arma::fill::zeros);
	for (int j = 1; j < size; j++) 
		output(j - 1) = log( y(j) / y(j - 1) ) / log( x(j) / x(j - 1) );
	return output;
}

template<typename _type>
_type simpson_rule(const arma::vec& x, const arma::Col<_type>& f) {
	const int N = (int)f.size() - 1;
	arma::vec h(N);
	for (int i = 0; i < N; i++)
		h(i) = x(i + 1) - x(i);
	
	_type sum = _type(0.0);
//#pragma omp parallel for reduction(+: sum)
	for (int i = 0; i <= N / 2 - 1; i++) {
		_type a = 2 - h(2 * i + 1) / h(2 * i);
		_type b = (h(2 * i) + h(2 * i + 1)) * (h(2 * i) + h(2 * i + 1)) / (h(2 * i) * h(2 * i + 1));
		_type c = 2 - h(2 * i) / h(2 * i + 1);
		sum += (h(2 * i) + h(2 * i + 1)) / 6.0 * (a * f(2 * i) + b * f(2 * i + 1) + c * f(2 * i + 2));
	}

	if (N % 2 == 0) {
		_type a = (2 * h(N - 1) * h(N - 1)  + 3 * h(N - 1) * h(N - 2)) / (6 * (h(N - 2) + h(N - 1)));
		_type b = (	h(N - 1) * h(N - 1) 	+ 3 * h(N - 1) * h(N - 2)) / (6 *  h(N - 2));
		_type c = (	h(N - 1) * h(N - 1) * h(N - 1)				 	 ) / (6 *  h(N - 2) * (h(N - 2) + h(N - 1)));
	}
	return sum;
}
DISABLE_WARNING_POP
template cpx simpson_rule<cpx>(const arma::vec&, const arma::Col<cpx>&);
template double simpson_rule<double>(const arma::vec&, const arma::Col<double>&);


/// <summary>
/// find non-unique elements in input array and store only elemetns, which did not have duplicates
/// </summary>
/// <param name="arr_in"> degenerated input array </param>
/// <returns> indices to unique elements </returns>
arma::vec get_NonDegenerated_Elements(const arma::vec& arr_in) {
	std::vector<double> arr_degen;
	u64 N = arr_in.size();
	NO_OVERFLOW(
		for (int k = 0; k < N - 1; k++)
			if (abs(arr_in(k + 1) - arr_in(k)) < 1e-12)
				arr_degen.push_back(arr_in(k));
	);
	if (abs(arr_in(N - 1) - arr_in(N - 2)) < 1e-12)
		arr_degen.push_back(arr_in(N - 1));

	std::vector<double> arr_unique = arma::conv_to<std::vector<double>>::from(arr_in);
	for (auto& E : arr_degen) {
		auto new_end = remove_if(arr_unique.begin(), arr_unique.end(), [&](double a) {
			return abs(a - E) < 1e-12;
			});
		arr_unique.erase(new_end, arr_unique.end());
	}
	return arma::unique((arma::vec)arr_unique);
}

/// <summary>
///
/// </summary>
/// <param name="arr_in"></param>
/// <returns></returns>
double binder_cumulant(const arma::vec& arr_in) {
	double av2 = 0, av4 = 0;
#pragma omp parallel for reduction(+: av2, av4)
	for (int k = 0; k < arr_in.size(); k++) {
		double val = arr_in(k) * arr_in(k);
		av2 += val;
		av4 += val * val;
	}
	return 1 - av4 * (double)arr_in.size() / (3.0 * av2 * av2);
}