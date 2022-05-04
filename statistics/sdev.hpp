#ifndef statistics_sdev_H
#define statistics_sdev_H

#include <vector>
#include <cmath>
#include <statistics/svar.hpp>
using namespace std;

template<class d>
d sdev(vector<d> &v)
{return sqrt(svar(v));}
#endif // statistics_sdev_H
	
