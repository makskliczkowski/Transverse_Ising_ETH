#pragma once
#include <random>
#include <numeric>
#include "xoshiro_pp.h"

#ifndef RANDOM_H
#define RANDOM_H
class random_num {
    // Hold RNG state as a member variable
    //std::mt19937 eng;
    XoshiroCpp::Xoshiro256PlusPlus eng;
    std::uniform_real_distribution<double> uni_double_dist;
    std::normal_distribution<double> normal_dist;
public:
    random_num(const std::uint64_t seed = static_cast<uint64_t>(std::time(nullptr))){

        //std::random_device rd{};
        //std::seed_seq seq{rd(), rd(),rd(),rd()};
        //this->eng = std::mt19937{ seq};//std::chrono::high_resolution_clock::now().time_since_epoch().count()} };
        this->eng = XoshiroCpp::Xoshiro256PlusPlus(this->Random_SeedInit(seed));
        //eng.discard(700000);

    }
    uint64_t Random_SeedInit(uint64_t n) const
    {
        std::vector<uint64_t> s(16, 0);
        for (int i = 0; i < 16; i++)
        {
            n ^= n >> 12;   // a
            n ^= n << 25;   // b
            n ^= n >> 27;   // c

            // 2685821657736338717 = 72821711 * 36882155347, from Pierre L'Ecuyer's paper
            s[i] = n * 2685821657736338717LL;
        }
        return std::accumulate(s.begin(), s.end(), 0.0);
    }

    inline double rand_uni_dist(double _min, double _max){
        //return uni_double_dist(_min, _max)(this->eng);
        //return XoshiroCpp::DoubleFromBits(this->eng());
        return std::uniform_real_distribution<double>(_min, _max)(this->eng);
    }
    inline int rand_uni_int_dist(int L){
        return std::uniform_int_distribution<>(0, L)(eng);
    }
    
    double rand_norm_dist() {
        return normal_dist(eng);
    }
};


#endif