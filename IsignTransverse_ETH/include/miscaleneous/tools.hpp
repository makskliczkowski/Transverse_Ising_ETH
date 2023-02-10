#pragma once
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <numeric>

// Calculate elapsed time from start
//<! start --> starting time_point to substract from current time and get distance
inline double tim_s(clk::time_point start) {
	return double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(\
		clk::now() - start)).count()) / 1000.0;
}
//-------------------------------------------------------------------------------------------------------------- OPERATION ON STRINGS
// checking if string is a number
inline bool isNumber(const std::string& str) {
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

// split string into pieces divided by 'delimeter'
// as output is vector of strings
inline std::vector<std::string> 
split_str(
	std::string s,			//<! string to split
	std::string delimiter	//<! symbol serving as split value
) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(s.substr(pos_start));
	return res;
}


// Print value of any type to custom precision
template <typename _type> 
inline std::string to_string_prec(
	const _type a_value, //<! argument to convert parse to stream
	const int n = 3		 //<! precision/number of digits after comma
) {
	std::ostringstream outie;
	outie.precision(n);
	outie << std::fixed << a_value;
	return outie.str();
}

//-------------------------------------------------------------------------------------------------------------- PERMUTATION AND SORTING OF DATA

// Sorts the vector and saves the permutation with a lambda like function compare
template <has_access_operator container, callable_type F> 
inline std::vector<std::size_t> sort_permutation(
	const container& vec,	//<! input vector to find permutation on
	F&& compare						//<! callable/invocable
) {
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) {
			return compare(vec[i], vec[j]);
		});
	return p;
}

// Applies permutation on a given vector
template <has_access_operator container> 
inline void apply_permutation(
	container& vec,			//<! vector to permute
	const std::vector<std::size_t>& p	//<! permutation on input vector
) {
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

//-------------------------------------------------------------------------------------------------------------- ADDITIONAL TOOLS
template <typename T> 
inline int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}


// Overriding the ostream operator for pretty printing vectors.
template <has_output_operator _type> // ensure output stream for object inside container
inline std::ostream& operator<<(
	std::ostream& os,				//<! output stream
	const std::vector<_type>& vec	//<! vector to print
) {
	if (vec.size() != 0) {
		std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<_type>(os, " "));
		os << vec.back() << ' ';
	}
	else
		os << "Empty container!" << std::endl;
	return os;
}