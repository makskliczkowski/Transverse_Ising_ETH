#pragma once
#ifndef RANDOM_H
#define RANDOM_H

#include "../src/xoshiro_pp.h"
#include <random>
#include <ctime>
#include <numeric>

/// <summary>
/// Random number generator class
/// </summary>
class randomGen {
private:
	XoshiroCpp::Xoshiro256PlusPlus engine;
	std::uint64_t init_seed;
public:
	explicit randomGen(const std::uint64_t seed = std::random_device{}()) {
		this->init_seed = seed;
		this->engine = XoshiroCpp::Xoshiro256PlusPlus(this->SeedInit(seed));
	}
	uint64_t SeedInit(uint64_t n) const
	{
		std::vector<uint64_t> s(16, 0);
		for (int i = 0; i < 16; i++)
		{
			n ^= n >> 12;   // a
			n ^= n << 25;   // b
			n ^= n >> 27;   // c
			s[i] = n * 2685821657736338717LL;                               // 2685821657736338717 = 72821711 * 36882155347, from Pierre L'Ecuyer's paper
		}
		return std::accumulate(s.begin(), s.end(), 0.0);
	}
	void reset()
		{this->engine = XoshiroCpp::Xoshiro256PlusPlus(this->SeedInit(this->init_seed));}

	//------------------------------------------------------------------------------ WRAPPERS ON RANDOM FUNCTIONS
	std::complex<double> randomCpx_uni(double _min = 0, double _max = 1) 
	{
		std::uniform_real_distribution<double> dist(_min, _max);
		return std::complex<double>(dist(engine), dist(engine));
	}
	double randomReal_uni(double _min = 0, double _max = 1) 
		{ return std::uniform_real_distribution<double>(_min, _max)(engine); }
	uint64_t randomInt_uni(int _min, int _max) 
		{ return std::uniform_int_distribution<uint64_t>(_min, _max)(engine); }


	template <typename _type> 
	_type random_uni(double _min, double _max) 
		{ return std::uniform_real_distribution<_type>(_min, _max)(engine); }

	template <typename _type> 
	arma::Col<_type> create_random_vec(const uint64_t size, double h = 0.0) 
	{
		arma::Col<_type> random_vec(size, arma::fill::zeros);
		for (u64 j = 0; j <= size / 2.; j++) {
			u64 idx = size / (long)2 - j;
			random_vec(idx) = random_uni<_type>(-h, h);
			idx += 2 * j;
			if (idx < size) random_vec(idx) = random_uni<_type>(-h, h);
		}
		return arma::normalise(random_vec);
	}
	template <typename _type> 
	std::vector<_type> create_random_stdvec(const uint64_t size, double h = 0.0) 
	{
		std::vector<_type> random_vec(size, 0);
		_type norm = 0;
		for (uint64_t j = 0; j < size; j++) {
			random_vec[j] = random_uni<_type>(-h, h);
			norm += random_vec[j];
		}
		for (uint64_t j = 0; j < size; j++)
			random_vec[j] /= norm;
		return random_vec;
	}
};

template <> 
inline 
std::complex<double> 
	randomGen::random_uni<std::complex<double>>(
		double _min, double _max)
	{ return randomCpx_uni(_min, _max); }

	
#endif // !RANDOM_H
