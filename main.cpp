#include <statistics.hpp>
#include <iostream>
#include <vector>

using namespace std;

int main()
{
	vector<double> v1 = {1,2,3,4};
	cout << Mean(v1) << endl;
	cout << svar(v1) << endl;
	cout << sdev(v1) << endl;
}
