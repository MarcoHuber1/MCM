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

const char* Datei = "data1.txt";

double (*f)(double);
double (*p)(double, double, double);
f = Cos;
p = Gaussdist;

double a = integrate(p,f);
cout << "Solution of Integral with trap. approx:" << a << endl;

///////////////////////////////////



vector<double> N(101,0);
FILE * handle = fopen(Datei, "w");
for(int j = 0; j<101 ; ++j)
{
	
 vector<double> M(100,1);
 for(int m = 0; m<100; ++m)
 {
	N[j] = pow(10,1.+j*5/100.);
	fprintf(handle, "%ld ",lround(N[j]));
	//cout<< lround(N[j]) << endl;
	vector<double> v(lround(N[j]),0);
	RNG_normv(v);
	for(int i = 0; i < lround(N[j]); ++i)
	{
		v[i] = cos(v[i]);
	}
	M[m] = Mean(v);

	fprintf(handle, "%lf ",M[m]);
 }
	fprintf(handle, "\n");

}
fclose(handle);
}







