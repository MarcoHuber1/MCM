#ifndef autocorrelation_autocorrelationfuction_H
#define autocorrelation_autocorrelationfuction_H

template<typename Val>
double autocorrfunc(int t, std::vector<Val> &Y)
{
	int Dim = Y.size();
	double prefactor = (double)1/(Dim-t);
	double rho_t = 0;
	double C_t = 0;
	double C_0 = 0;
	
	double y_minus = 0;
	double y_plus = 0;
	
	// calculation of y+ and y-
	for(int i = 0; i<(Dim-t); ++i)
	{
		y_minus += Y[i];
		y_plus += Y[i+t];
	}
	y_minus *= prefactor;
	y_plus *= prefactor;

	//std::cout << t << " " << y_minus << " " << y_plus << std::endl;
	
	//calculation of C_0
	for(int j = 0; j<Dim; ++j)
	{C_0 += (Y[j] - y_minus) * (Y[j] - y_plus); }
	C_0 *= (double)1/Dim;
	
	//calculation of C_t
	for(int j = 0; j<(Dim-t); ++j)
	{C_t += (Y[j] - y_minus) * (Y[j+t] - y_plus); }
	C_t *= prefactor;
	
	rho_t = C_t/C_0;
	//std::cout << rho_t <<std::endl;
	return rho_t;
	
}
		
		
template<typename Val>
double tau_int(std::vector<Val> &V)
{
	double tau = 0.5;
	int t = 0;
	while(true)
	{
		double Contribution = autocorrfunc(t, V);
		if(Contribution < 0)
		break;
		
		tau += Contribution;
		t += 1;
		//std:: cout << t << " " <<  Contribution <<std::endl;
	}
	std::cout << tau << std::endl;
	return tau;
}

template<typename Val>
double auto_std_err_prim(int Dim, double tau, std::vector<Val> &V)
{
    return sqrt(2*tau/Dim)*sdev_s(V, Mean(V));
}

#endif // autocorrelation_autocorrelationfuction_H
