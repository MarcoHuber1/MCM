#include <iostream>
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

template<class d>
auto sdev(vector<d> &v)
{
	d meanvalue = Mean(v);
}


int main()
{
vector<double> v = {1,2.3,3.4,3};
cout << Mean(v) << endl;
}
