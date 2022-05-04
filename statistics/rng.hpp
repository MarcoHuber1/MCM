#ifndef statistics_rng_H
#define statistics_rng_H

// mersenne_twister_engine:: seed, uniform distribution, normal distribution
#include <iostream>
#include <chrono>
#include <random>
#include <vector>

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 gen;

//uniform numbers

template<typename d>
void RNG_univ(vector<d> &v, unsigned sd = seed)
{
	gen.seed(sd);
	std::uniform_int_distribution<int> unidist(1,100);

        for(int i=0; i<v.size(); ++i)
        {
                v[i] = unidist(gen);
        }
}

auto RNG_uni(unsigned sd = seed)
{
	gen.seed(sd);
	std::uniform_int_distribution<int> unidist(1,100);
	return unidist(gen);
}


//normal numbers

template<typename d>
void RNG_normv(vector<d> &v, double mean = 0.0, double sdev = 1.0, unsigned sd = seed)
{
        gen.seed(sd);
        std::normal_distribution<double> normdist(mean,sdev);

        for(int i=0; i<v.size(); ++i)
        {
                v[i] = normdist(gen);
        }
}

auto RNG_norm(double mean = 0.0, double sdev = 1.0, unsigned sd = seed)
{
        gen.seed(sd);
        std::normal_distribution<double> normdist(mean,sdev);
        return normdist(gen);
}


#endif // statistics_rng_H

