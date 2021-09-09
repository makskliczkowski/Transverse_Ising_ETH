#include "include/headers.h"
int num_of_threads = 8;

std::random_device rd;
std::mt19937::result_type seed = 87178291199L; // set constant to maintain same disorder for different sizes etc
// rd() ^ (\
	(std::mt19937::result_type)\
	std::chrono::duration_cast<std::chrono::seconds>(\
		std::chrono::system_clock::now().time_since_epoch()\
		).count() +\
	(std::mt19937::result_type)\
	std::chrono::duration_cast<std::chrono::microseconds>(\
		std::chrono::high_resolution_clock::now().time_since_epoch()\
		).count());
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
double simpson_rule(double a, double b,
	int n, // Number of intervals
	const arma::vec& f)
{
	double h = (b - a) / n;

	// Internal sample points, there should be n - 1 of them
	double sum_odds = 0.0;
#pragma omp parallel for reduction(+: sum_odds)
	for (int i = 1; i < n; i += 2) {
		int idx = ((a + i * h) + abs(a)) / h;
		sum_odds += f(idx);
	}
	
	double sum_evens = 0.0;
#pragma omp parallel for reduction(+: sum_evens)
	for (int i = 2; i < n; i += 2) {
		int idx = ((a+i * h) + abs(a)) / h;
		sum_evens += f(idx);
	}
	
	return (f(0) + f(f.size() - 1) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}