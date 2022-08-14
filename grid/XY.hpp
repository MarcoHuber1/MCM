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
        vec[point] = 0;
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
        for(int neighbor = 0; neighbor<4; ++neighbor)
        {
            Energy -= cos(Theta[point])*cos(Theta[g->getNN(point,neighbor)]) + sin(Theta[point])*sin(Theta[g->getNN(point,neighbor)]);
        }
    }
    return 0.5*g->getJ()*Energy/g->Dim();
}
//Magnetisationdensity
template<typename Val>
double MD_XY(Spin_vector<Val> &Theta, Grid_XY *g) //Energy density
{
    double Magnetization = 0;

    for(int point = 0; point< Theta.Dim(); ++point)
    {
        Magnetization += Theta[point];
        
    }
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
void Guidance_Calc_i(std::mt19937 gen, Spin_vector<Val> p, Spin_vector<Val> Theta, double &p_con_squared, double &H_g, Grid_XY *g)
{
    //Conjugate Momenta and Guidance Hamiltonian
    std::normal_distribution<> nd{0,1};
    
    for(int i = 0; i < Theta.Dim(); ++i)
    {p[i] = nd(gen);}

    for(int i = 0; i < Theta.Dim(); ++i)
    {p_con_squared += pow(p[i],2);}
    
    //Guidance Hamiltonian
    H_g = p_con_squared/2 + g->getBeta()*ED_XY(Theta,g);

}

template<typename Val>
void Guidance_Calc_f(Spin_vector<Val> Theta, double &p_con_squared, double &H_g, Grid_XY *g)
{
    //Guidance Hamiltonian
    H_g = p_con_squared/2 + g->getBeta()*ED_XY(Theta,g);
}

template<typename Val>
double dVdq(int i,Spin_vector<Val> Theta, Grid_XY *g)
{
    double dV = 0;

        for(int j = 0; j<4; ++j)
        {
            dV += sin(Theta[i]-Theta[g->getNN(i,j)]);

        }
    return g->getJ()* g->getBeta() * 0.5* dV;

}

template<typename Val>
void Leapfrog(Spin_vector<Val> &Theta, Spin_vector<Val> &p, int &t_LF, Grid_XY *g)
{
    double stepsize = 0.1;
    //q_i(0) = Configuration
    Spin_vector<Val> q(g);
    for(int m =0; m<Theta.Dim(); ++m)
    {
        q[m] = Theta[m];
    }

    //initial Halfstep
    for(int i = 0; i<p.Dim(); ++i)
    {
        p[i] -=dVdq(i,Theta,g) * stepsize/2;

    }

    //middle part
    for(int time = 0; time<(t_LF -1)*stepsize; ++time)
    {
        for(int i = 0; i<p.Dim(); ++i)
        {
            q[i] += p[i]*stepsize;
            p[i] -= dVdq(i,Theta,g)*stepsize;

        }
    }

    //final half step
    for(int i = 0; i<p.Dim(); ++i)
        {
            q[i] += p[i]*stepsize;

            Theta[i] = q[i];
            p[i] -= dVdq(i,Theta,g)*stepsize/2;
        }

}

template<typename Val>
void HMC(Grid_XY *g, Spin_vector<Val> &Theta, std::mt19937 &gen, int &t_F_Max, int &t_LF)
{
    for(int t_F = 0; t_F < t_F_Max; ++t_F)
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

        Leapfrog(Final, p_i_initial, t_LF,g);

        Guidance_Calc_f(Final,p_con_squared_initial,H_g_final,g);

        //Accept reject method
        //double dH
        double acceptance = exp(H_g_initial - H_g_final);
        double RN = 0;
        RNG_uni(RN,gen);

        //accept new config
        if(RN <= acceptance)
        {
            for(int index = 0; index < Theta.Dim(); ++index)
            {Theta[index] = Final[index];}
        }

    }
    
}

#endif //Grid_XY_H
