#include <iostream>
#include <statistics.hpp>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES

using namespace std;

//trapezoidal integral
double integrate(double (*p)(double, double, double),double (*f) (double),double a = -1000, double b = 1000, int steps = 10000)
{	
	double result = 0.5 * (f(a)*p(a,0,1) + f(b)*p(b,0,1));
	double stepsize = (b - a)/steps;

	for(int i=1; i < steps; ++i)
	{result += (f(a + i*stepsize) * p(a + i*stepsize,0,1));}
	result *= stepsize;
	return result;
}

double Gaussdist(double x, double mean = 0, double sdev = 1)
{return (1/(sdev * sqrt(2 * M_PI)) * exp(-0.5 * pow((x-mean)/(sdev),2)));}

auto Cos(double x)
{return cos(x);}

int main()
{

double (*f)(double);
double (*p)(double, double, double);
f = Cos;
p = Gaussdist;

double a = integrate(p,f);
double b = Gaussdist(1);
cout << b << endl;
cout << a << endl;

}
