#include <iostream>
#include <statistics.hpp>
#include <vector>
#include <cmath>
#include <iomanip>

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

const char* Datei = "IsingE.txt";
const char* Datei2 = "IsingM.txt";

double (*f)(double);
double (*p)(double, double, double);
f = Cos;
p = Gaussdist;

double a = integrate(p,f);
cout << "Solution of Integral with trap. approx:" << a << endl;

///////////////////////////////////
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 gen;
gen.seed(seed);

//////////////////////////////////

vector<double> N(101,0);
FILE * handle = fopen(Datei, "w");
FILE * handle2 = fopen(Datei2, "w");

for(int j = 0; j<101 ; ++j) // 101 different Ns
{ 
 N[j] = pow(10,1.+j*5/100.); //sample size N

 fprintf(handle, "%ld ",lround(N[j]));
 fprintf(handle2, "%ld ",lround(N[j]));

 vector<double> M(100,0); //vector for Means of rep. sampling
 vector<double> stdev_vector(101,0); //standard dev. values

 for(int m = 0; m<100; ++m) //repeated sampling
 {
	vector<double> v(lround(N[j]),0);
	RNG_normv(v, gen);
	for(int i = 0; i < lround(N[j]); ++i) //Samples i=1,...,N
	{
		v[i] = cos(v[i]); 
	}
	M[m] = Mean(v);

	fprintf(handle, "%lf ",M[m]);
	
	stdev_vector[j] = sdev_s(M,a);
 }

	fprintf(handle2, "%lf ",stdev_vector[j]);
	fprintf(handle, "\n");
	fprintf(handle2, "\n");

}
fclose(handle);
fclose(handle2);
}







