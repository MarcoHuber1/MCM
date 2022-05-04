#include "mean.hpp"

template<class d>
auto Mean(vector<d> &v)
{
	d meanvalue = 0;
	for(auto i:v)
	{
		meanvalue += i;
	}
	meanvalue = meanvalue/v.size();
	return meanvalue;
}


