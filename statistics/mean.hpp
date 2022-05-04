#ifndef statistics_Mean_H
#define statistics_Mean_H
#include <vector>
using namespace std;

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
#endif // statistics_Mean_H


