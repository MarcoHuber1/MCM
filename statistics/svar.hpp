#ifndef statistics_svar_H
#define statistics_svar_H
#include <vector>
#include <cmath>
#include <statistics/mean.hpp>
using namespace std;

template<class d>
auto svar(vector<d> &v)
{
	d meanvalue = Mean(v);
        d var = 0;

        for(auto i:v)
        {
                var += pow((i - meanvalue),2);
        }
        var = var/v.size();
        return var;
}

template<class d>
auto svar_s(vector<d> &v,d meanvalue)
{
        d var = 0;

        for(auto i:v)
        {
                var += pow((i - meanvalue),2);
        }
        var = var/(v.size()-1);
        return var;
}


#endif // statistics_svar_H

