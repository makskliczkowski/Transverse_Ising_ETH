#include "include/headers.h"


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
    for (char const &c : str) {
		if(c == '.'&& !found_dot ){
			found_dot = true;
			continue;
		}
		if(c=='-' && !found_minus){
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

