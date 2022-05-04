#include <vector>
#include <cmath>

using namespace std;

template<class d>
auto svar(vector<d> &v)
{
        d meanvalue = Mean(v);
        d dev = 0;
        for(auto i:v)
        {
                dev += pow((i - meanvalue),2);
        }
        dev = dev/v.size();
        return dev;
}

