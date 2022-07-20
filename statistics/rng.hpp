#ifndef statistics_rng_H
#define statistics_rng_H

// mersenne_twister_engine:: seed, uniform distribution, normal distribution
#include <iostream>
#include <chrono>
#include <random>
#include <vector>

//uniform numbers

template<typename d>
void RNG_univ(vector<d> &v, std::mt19937 &gen)
{
	std::uniform_int_distribution<int> unidist(1,100);

        for(int i=0; i<v.size(); ++i)
        {
                v[i] = unidist(gen);
        }
}
template<typename d>
void RNG_uni(d &number, std::mt19937 &gen)
{
	std::uniform_real_distribution<> unidist(0,1);
	number = unidist(gen);
}


//normal numbers

template<typename d>
void RNG_normv(vector<d> &v, std::mt19937 &gen, double mean = 0.0, double sdev = 1.0)
{
        std::normal_distribution<double> normdist(mean,sdev);

        for(int i=0; i<v.size(); ++i)
        {
                v[i] = normdist(gen);
        }
}

template<typename d>
void RNG_norm(d &number, std::mt19937 &gen, double mean = 0.0, double sdev = 1.0)
{
        std::normal_distribution<d> normdist(mean,sdev);
        number = normdist(gen);
}


#endif // statistics_rng_H

