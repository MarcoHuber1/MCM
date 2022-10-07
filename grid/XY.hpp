#ifndef Grid_XY_H
#define Grid_XY_H
#include<statistics.hpp>
#include<autocorrelation.hpp>
#include<cstddef>
#include<array>
#include<math.h>


/////////////////////////////Grid_XYClass/////////////////////////
class Grid_XY
{
	private:
	size_t m_dim_x = 10;
	size_t m_dim_y = 10;
	size_t m_Dim = m_dim_x * m_dim_y;

    double mJ = 1;
    double mk = 1;
    double mB = 1;
    double mbeta = 1/mT;
    double mT = 1;

    std::vector<std::array<int,4>> NextNeigbor;


	public:	
	Grid_XY(size_t Lx, size_t Ly);
	Grid_XY(size_t L);


    //Functions
    //Get Grid_XY dimensions
    size_t Dim() {return m_Dim;};
    size_t DimX() {return m_dim_x;};
    size_t DimY() {return m_dim_y;};

    //parameters
    void setJ(double newJ) {mJ = newJ;};
    void setk(double newk) {mk = newk;};
    void setB(double newB) {mB = newB;};
    void setBeta(double newBeta) {mbeta = newBeta;};
    void setT(double newT) {mT = newT;};

    double getJ() {return mJ;};
    double getk() {return mk;};
    double getB() {return mB;};
    double getBeta() {return mbeta;};
    double getT() {return mT;};


    //Set NN
    void setNN(int i, int j, int value){NextNeigbor[i][j] = value;};
    //void resizeNN(){NextNeigbor.resize(m_Dim);}
    double getNN(int i, int j){return NextNeigbor[i][j];};

};

Grid_XY::Grid_XY(size_t L)
{
	m_dim_x = L;
	m_dim_y = L;
	m_Dim = m_dim_x * m_dim_y;
    NextNeigbor.resize(m_Dim);
}

Grid_XY::Grid_XY(size_t Lx, size_t Ly)
{
	m_dim_x = Lx;
	m_dim_y = Ly;
	m_Dim = m_dim_x * m_dim_y;
    NextNeigbor.resize(m_Dim);
}


/////////////////////////////GRIDFUNCTIONS/////////////////////////


/////////////////////////////Spin_vectorClass/////////////////////////
template<typename Val>
class Spin_vector //1D Grid_XY representation
{
    private:
    std::vector<Val> vec;

    Grid_XY *grid; //get access to private of grid

    public:
    //Constructors
    Spin_vector(Grid_XY *grid);
    Spin_vector(Grid_XY *grid, std::mt19937 &gen);
		
    //Destructor
    ~Spin_vector();
		
    //Operators
    Val& operator[](size_t index);
        
    //Functions
    size_t Dim() {return grid->Dim();}    
	
    void print();

};

//Constructors

template<typename Val>
Spin_vector<Val>::Spin_vector(Grid_XY *grid,std::mt19937 &gen) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);
    auto pi = M_PI;
    std::uniform_real_distribution<> unidist(0,2*pi);

     for(int point = 0; point< DIM; ++point)
    {
        vec[point] = unidist(gen);
    }

}
template<typename Val>
Spin_vector<Val>::Spin_vector(Grid_XY *grid) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);
    for(int point = 0; point< DIM; ++point)
    {
        vec[point] = 1.5*M_PI;
    }

}

//Destructor
template<typename Val>
Spin_vector<Val>::~Spin_vector()
{}

//operators
template<typename Val>
Val& Spin_vector<Val>::operator[](size_t index)
{
	return vec[index];
}

/////////////////////////////VECTORFUNCTIONS/////////////////////////
template<typename Val>
void Spin_vector<Val>::print() //prints whole vector
{
	auto D = Dim();

	for(int d = 0; d<D; ++d)
	{
		std::cout << vec[d] << " ";
	}
	std::cout << "\n";
    std::cout << "\n \n";
}
//Energydensity
template<typename Val>
double ED_XY(Spin_vector<Val> &Theta, Grid_XY *g) //Energy density
{
    double Energy = 0;

    for(int point = 0; point<Theta.Dim(); ++point)
    {
    /*
        for(int neighbor = 0; neighbor<4; ++neighbor)
        {
            Energy -= cos(Theta[point] - Theta[g->getNN(point,neighbor)]);
        }
        //(cos(Theta[point])*cos(Theta[g->getNN(point,neighbor)]) + sin(Theta[point])*sin(Theta[g->getNN(point,neighbor)]));
   */
        Energy -= cos(Theta[point] - Theta[g->getNN(point,1)]) + cos(Theta[point] - Theta[g->getNN(point,3)]);
    }
    return 0.5*g->getJ()*Energy/g->Dim();
}
//Magnetisationdensity
template<typename Val>
double MD_XY(Spin_vector<Val> &Theta, Grid_XY *g) //Energy density
{
    double Magnetization = 0;
    double spinx = 0;
    double spiny = 0;

    for(int point = 0; point< Theta.Dim(); ++point)
    {
        spinx += cos(Theta[point]);
        spiny += sin(Theta[point]);
    }
    Magnetization = sqrt(pow(spinx,2) + pow(spiny,2));
    return abs(Magnetization)/g->Dim();
}


/////////////////////////////Lattice_XYClass/////////////////////////
template<typename Val>
class Lattice_XY //2D representation of Grid_XY
{
    private:
    std::vector<std::vector<Val>> Lat;

    Grid_XY *grid;
    
    public:
    Lattice_XY(Grid_XY *grid);
    
    ~Lattice_XY();

    Val& operator()(size_t i, size_t j);

    size_t DimX() {return grid->DimX();}
    size_t DimY() {return grid->DimY();}
    
    void print();
};

//Constructors
template<typename Val>
Lattice_XY<Val>::Lattice_XY(Grid_XY *grid) :grid(grid)
{
    auto X = grid->DimX();
    auto Y = grid->DimY();
    
    Lat.resize(X);
    for(int i=0; i<X; ++i)
    {
        Lat[i].resize(Y);
    }
    
    for(int i = 0; i<X; ++i)
	{
		for(int j = 0; j<Y; ++j)
		{
			Lat[i][j] = i*Y+j;
		}
    }
}

//Destructors
template<typename Val>
Lattice_XY<Val>::~Lattice_XY()
{}

//Operators
template<typename Val>
Val& Lattice_XY<Val>::operator()(size_t x, size_t y)
{
	return Lat[x][y];
}
/////////////////////////////LATTICEFUNCTIONS/////////////////////////
template<typename Val>
void Lattice_XY<Val>::print()
{
	auto X = DimX();
	auto Y = DimY();

	for(int i = 0; i<X; ++i)
	{
		for(int j = 0; j<Y; ++j)
		{
			std::cout << Lat[i][j] << " ";
		}
	 	std::cout << "\n";
    }
	std::cout << "\n \n";
}


template<typename Val>
Spin_vector<Val> transform(Lattice_XY<Val> &lat, Grid_XY *grid)
{
	Spin_vector<Val> vec(grid);

	auto X = lat.DimX();
	auto Y = lat.DimY();

	for(int i = 0; i<X; ++i)
	{
		for(int j=0; j<Y; ++j)
		{
			vec[i*Y+j] = lat(i,j);
		}
	}

	return vec;
}

template<typename Val>
Lattice_XY<Val> transform(Spin_vector<Val> &vec, Grid_XY *grid) //transforms vec to lat
{
	Lattice_XY<Val> lat(grid);

	auto D = vec.Dim();
	auto X = grid->DimX();
	auto Y = grid->DimY();

         for(int i = 0; i<X; ++i)
         {
                 for(int j=0; j<Y; ++j)
                 {
                         lat(i,j) = vec[i*Y+j];
                 }
         }

	return lat;
}

/////////////////////////////////NextNeigbor/////////////////////
template<typename Val>
void NN(Grid_XY *g, Spin_vector<Val> &Configuration)
{
    //NextNeigbor = Spin_vector(top,bottom,left,right)

    Lattice_XY<Val> lat(g); //lat(x,y)

    int x,y = 0;

    const auto Lx = g->DimX();
    const auto Ly = g->DimY();

    for(int point = 0; point < g->Dim(); ++point)
    {
        x = point/Ly; //X position in Grid_XY (up to down)
        y = point % Ly;    //Y position in Grid_XY (left to right)

        for(int i = 0; i<4; ++i)
        {
            if(i == 0) //top
            {g->setNN(point,i,lat((x-1+Lx)%Lx,y));}

            if(i == 1)//bottom
            {g->setNN(point,i,lat((x+1)%Lx,y));}

            if(i == 2) //left
            {g->setNN(point,i,lat(x,(y+Ly-1)%Ly));}

            if(i == 3) //right
            {g->setNN(point,i,lat(x,(y+1)%Ly));}
        }
    }
}


//HMC-Method//////////////////////////////////////////////////
template<typename Val>
void Guidance_Calc_i(std::mt19937 gen, Spin_vector<Val> &p, Spin_vector<Val> &Theta, double &p_con_squared, double &H_g, Grid_XY *g)
{
    //Conjugate Momenta and Guidance Hamiltonian
    std::normal_distribution<> nd{0,1};

    for(int i = 0; i < Theta.Dim(); ++i)
    {p[i] = nd(gen);}

    for(int i = 0; i < Theta.Dim(); ++i)
    {p_con_squared += pow(p[i],2);}

    //Guidance Hamiltonian
    H_g = p_con_squared/2 + g->getBeta()*ED_XY(Theta,g)*g->Dim();

}

template<typename Val>
void Guidance_Calc_f(Spin_vector<Val> &Theta, Spin_vector<Val> &p, double &p_con_squared, double &H_g, Grid_XY *g)
{
    for(int i = 0; i < Theta.Dim(); ++i)
    {
        p_con_squared += pow(p[i],2);
    }

    //Guidance Hamiltonian
    H_g = p_con_squared/2 + g->getBeta()*ED_XY(Theta,g)*g->Dim();
}

template<typename Val>
double dVdq(int i,Spin_vector<Val> &Theta, Grid_XY *g)
{
    double dV = 0;

        for(int j = 0; j<4; ++j)
        {
            dV += sin(Theta[i]-Theta[g->getNN(i,j)]);

        }

    return g->getJ()* g->getBeta() *dV;

}

template<typename Val>
void Leapfrog(Spin_vector<Val> &Theta, Spin_vector<Val> &p, double &t_LF, Grid_XY *g)
{
    double stepsize = 0.01;

    Spin_vector<Val> q(g);
    double steps = t_LF/stepsize;


    for(int m =0; m<Theta.Dim(); ++m)
    {
        q[m] = Theta[m];
    }

    //initial Halfstep (deltat/2)
    for(int i = 0; i<p.Dim(); ++i)
    {
        p[i] -=dVdq(i,q,g) * stepsize/2;
    }

    //middle part (t+deltat/2)
    for(double time = 0; time<(steps-1)*stepsize; time+=stepsize)
    {
        for(int i = 0; i<q.Dim(); ++i)
        {
            q[i] += p[i]*stepsize;
        }

        for(int i = 0; i<p.Dim(); ++i)
        {
            p[i] -= dVdq(i,q,g)*stepsize;
        }


    }

    //final half step
    for(int i = 0; i<p.Dim(); ++i)
    {
        q[i] += p[i]*stepsize;
        Theta[i] = q[i];
    }

    for(int i = 0; i<p.Dim(); ++i)
    {
       p[i] -= dVdq(i,q,g)*stepsize/2;
    }


}

template<typename Val>
void HMC(Grid_XY *g, Spin_vector<Val> &Theta, std::mt19937 &gen, int &t_HMC, double &t_LF, double &n)
{
    std::vector<double> Energies;
    std::vector<double> Magnetizations;
    double EDavg = 0;
    double MDavg = 0;

    for(int t = 0; t < t_HMC; ++t)
    {

        //generalized Coord. q are the Values of Spin_vector Theta
        //Conjugate Momenta and Guidance Hamiltonian
        Spin_vector<Val> p_i_initial(g);
        Spin_vector<Val> p_i_final(g);

        double p_con_squared_initial = 0;
        double p_con_squared_final = 0;
        double H_g_initial = 0;
        double H_g_final = 0;

        //Initial values:
        Guidance_Calc_i(gen,p_i_initial,Theta,p_con_squared_initial,H_g_initial,g);


        //Make new vector for final proposal
        Spin_vector<Val> Final(g);
        for(int index = 0; index < Theta.Dim(); ++index)
        {Final[index] = Theta[index];}

        //Leapfrog Algo
        //std::cout << Final[0] << " " << Final[1] << std::endl;
        Leapfrog(Final, p_i_initial, t_LF,g);
        //std::cout << Final[0] << " " << Final[1] << std::endl;
        //std::cout << "LOL" << std::endl;

        Guidance_Calc_f(Final,p_i_initial,p_con_squared_final,H_g_final,g);


        //std::cout << "Guidance: " << H_g_initial << " " << H_g_final << std::endl;


        //Accept reject method
        double acceptance = std::min(1.0,exp(H_g_initial - H_g_final));

        if(H_g_final >= H_g_initial)
        {
            for(int index = 0; index < Theta.Dim(); ++index)
            {
                Theta[index] = Final[index];
                //std::cout << "accepted: " << acceptance <<std::endl;
            }
        }
        else
        {
            double RN = 0;
            RNG_uni(RN,gen);
            //std::cout << acceptance << std::endl;
            //accept new config
            if(RN <= acceptance)
            {
                for(int index = 0; index < Theta.Dim(); ++index)
                {
                    Theta[index] = Final[index];
                    //std::cout << "accepted: " << RN << " " << acceptance <<std::endl;
                }
            }
        }

        //Measurements:
        if(t > t_HMC/2)
        {
            //std::cout << ED_XY(Theta,g) << std::endl;
            Energies.push_back(ED_XY(Theta,g));
            Magnetizations.push_back(MD_XY(Theta,g));
        }

    }
    EDavg = Mean(Energies);
    MDavg = Mean(Magnetizations);
    //double tau_M_XY = tau_int(Magnetizations);
    //std::cout << g->getT() << " " << 2*tau_M_XY*15 << std::endl;
    //std::cout << n/0.1 << " " << 2*tau_M_XY*(n/0.1) << std::endl;
    std::cout <<  g->getT() << "  "<<EDavg << " " <<MDavg << std::endl;

    
}








//Delta Energy
template<typename Val>
void q_XY(double &q, Spin_vector<Val> &Theta, Grid_XY *g, int &point)
{
    for(int i = 0; i<4; ++i)
    {q -= cos(Theta[point] - Theta[g->getNN(point,i)]);}

}

//Spinflip
template<typename Val>
void Spinflip_proposal(Spin_vector<Val> &Theta,int &point, double &dTheta)
{
    Theta[point] += dTheta;
}

template<typename Val>
void Spinflip_back(Spin_vector<Val> &Theta,int &point, double &dTheta)
{
    Theta[point] -=dTheta;
}

template<typename Val>
void Metropolis_XY(Spin_vector<Val> &Theta, Grid_XY *g, int Iterations, std::mt19937 &gen)
{
    int MarkovTime_XY = 0;
/*
    const char* Datei = "IsingE100.txt";
    const char* Datei2 = "IsingM100.txt";
    FILE * handle = fopen(Datei, "w");
    FILE * handle2 = fopen(Datei2, "w");
*/

    double exponential_XY = 0;

    double VarianceE_XY = 0;
    double VarianceM_XY = 0;

    vector<double> ED_vector_XY;
    vector<double> MD_vector_XY;

    vector<double> ED_vector_squared_XY;
    vector<double> MD_vector_squared_XY;


    double MDavg_XY = 0;
    double EDavg_XY = 0;

    //Markov Process
    for(int i=0; i< Iterations; ++i)
    {
         //std::cout << "lol" << std::endl;
        //Getting each point of Configuration
        for(int point = 0; point < Theta.Dim(); ++point)
        {
            double Q_i  = 0; //Energy before
            double Q_f  = 0; //Energy after
            double dE = 0; //Energy difference


            //calculate Energy before spinflip
            q_XY(Q_i,Theta, g, point);

            //make spinflip with normal distr.
            double delta = 0.8;
            double dTheta = 0;
            RNG_norm(dTheta,gen,0,delta);
            Spinflip_proposal(Theta,point,dTheta); //new proposal Theta = Theta + dTheta

            //calculate Energy after spinflip
            q_XY(Q_f,Theta, g, point);

            dE = Q_f - Q_i;
            //std::cout << dE << std::endl;

            if(dE>0)
            {
                double RandomNumber = 0;
                RNG_uni(RandomNumber,gen);

                double Rho = exp(g->getBeta()*-dE);

                if(RandomNumber > Rho)
                {Spinflip_back(Theta,point,dTheta);}

            }
        }

        //fprintf(handle, "%lf\n",ED(Configuration,g));
        //fprintf(handle2, "%lf\n",MD(Configuration,g));
        MarkovTime_XY += 1;

        if(MarkovTime_XY > 5000) //Expectationvalue of E and M after equilibration and after each Markov Chain
        {
            ED_vector_XY.push_back(ED_XY(Theta,g));
            MD_vector_XY.push_back(MD_XY(Theta,g));
            //std::cout << "lol" << std::endl;

            ED_vector_squared_XY.push_back(pow(ED_XY(Theta,g),2));
            MD_vector_squared_XY.push_back(pow(MD_XY(Theta,g),2));
        }
    }
    //std::cout << "lal" << std::endl;

    MDavg_XY = Mean(MD_vector_XY);
    EDavg_XY = Mean(ED_vector_XY);
/*
    vector<double> C_eff_XY;
    vector<double> X_eff_XY;

    for(int m = 0; m < ED_vector_XY.size(); ++m)
    {
        C_eff_XY.push_back(pow(g->getBeta(),2)* Theta.Dim() * (ED_vector_squared_XY[m] - 2 * EDavg_XY * ED_vector_XY[m]));
        X_eff_XY.push_back(g->getBeta() * Theta.Dim() * (MD_vector_squared_XY[m] - 2 * MDavg_XY * MD_vector_XY[m]));
    }


    double tau_E_XY = tau_int(ED_vector_XY);
    double tau_E_squared_XY = tau_int(ED_vector_squared_XY);

    double tau_M_XY = tau_int(MD_vector_XY);
    double tau_M_squared_XY = tau_int(MD_vector_squared_XY);

    double tau_C_XY = tau_int(C_eff_XY);
    double tau_X_XY = tau_int(X_eff_XY);

    double sigma_E_XY = auto_std_err_prim(ED_vector_XY.size(),tau_E_XY, ED_vector_XY);
    double sigma_M_XY = auto_std_err_prim(MD_vector_XY.size(),tau_M_XY, MD_vector_XY);
    double sigma_C_XY = auto_std_err_prim(C_eff_XY.size(),tau_C_XY, C_eff_XY);
    double sigma_X_XY = auto_std_err_prim(X_eff_XY.size(),tau_X_XY, X_eff_XY);

    double C_XY = pow(g->getBeta(),2)* Theta.Dim() * svar(ED_vector_XY);
    double X_XY = g->getBeta() * Theta.Dim() * svar(MD_vector_XY);

    //std::cout << tau_C << " " << tau_X << std::endl;
    //std::cout << "tau C: "<<auto_std_err_prim(C_eff.size(),tau_C, C_eff) << std::endl;
    //std::cout << "tau X: "<<auto_std_err_prim(X_eff.size(),tau_X, X_eff) << std::endl;
    //std::cout << "tau MD: "<< auto_std_err_prim(MD_vector.size(),tau_M, MD_vector) << std::endl;
    //std::cout << "tau MD: "<< auto_std_err_prim(ED_vector.size(),tau_M, ED_vector) << std::endl;
    //std::cout << Blocking(MD_vector,20, g->getBeta() * Configuration.Dim()) << std::endl;
    //std::cout << Blocking(MD_vector,20, g->getBeta() * Configuration.Dim()) << std::endl;

    //std::cout << Bootstrap(MD_vector,tau_M, 1000, gen, g->getBeta()  * Configuration.Dim()) << std::endl;
    //std::cout << Bootstrap(MD_vector,tau_X, 1000, gen, g->getBeta()  * Configuration.Dim()) << std::endl;
    //std::cout << g->getT() << " " << EDavg << " " << MDavg << " " << C <<  " " << X << std::endl;
    //std::cout << g->getT() << " " <<sigma_E << " " << sigma_M << " " << sigma_C <<  " " << sigma_X << std::endl;
    //std::cout << g->getT() << " " <<tau_E << " " << tau_M << " " << tau_C <<  " " << tau_X << std::endl;
*/
    std::cout << g->getT() << " " <<EDavg_XY << " " << MDavg_XY <<std::endl;


}

#endif //Grid_XY_H
