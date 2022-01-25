#pragma once
#ifndef RANDOM_H
#define RANDOM_H

//#include "xoshiro_pp.h"
#include <random>
#include <ctime>
#include <numeric>

/// <summary>
/// Random number generator class
/// </summary>
class randomGen {
private:
	//XoshiroCpp::Xoshiro256PlusPlus engine; 
	std::mt19937_64 engine;
	std::mt19937::result_type seed = std::mt19937::result_type(87178291199L);
public:
	randomGen(std::uint64_t seed = std::random_device{}()) {
		//this->engine = XoshiroCpp::Xoshiro256PlusPlus(this->SeedInit(seed));
		this->engine = std::mt19937_64(seed);
	}
	uint64_t SeedInit(uint64_t n) const
	{
		uint64_t z = (n += 0x9e3779b97f4a7c15);
		z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
		z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
		return z ^ (z >> 31);
	}
	void reset() {
		this->engine = std::mt19937_64(seed);
	}

	/* WRAPPERS ON RANDOM FUNCTIONS */
	double randomReal_uni(double _min = 0, double _max = 1) {
		return std::uniform_real_distribution<double>(_min, _max)(engine);
	}
	uint64_t randomInt_uni(int _min, int _max) {
		return std::uniform_int_distribution<uint64_t>(_min, _max)(engine);
	}
};

#endif // !RANDOM_H