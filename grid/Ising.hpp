#ifndef Grid_H
#define Grid_H
#include<statistics.hpp>
#include<autocorrelation.hpp>
#include<cstddef>
#include<array>


/////////////////////////////GridClass/////////////////////////
class Grid
{
	private:
	size_t m_dim_x = 10;
	size_t m_dim_y = 10;
	size_t m_Dim = m_dim_x * m_dim_y;

    double mJ = 1;
    double mk = 1;
    double mB = 0;
    double mT = 1;
    double mbeta = 1/mT;


    std::vector<std::array<int,4>> NextNeigbor;


	public:	
	Grid(size_t Lx, size_t Ly);
	Grid(size_t L);


    //Functions
    //Get Grid dimensions
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

Grid::Grid(size_t L)
{
	m_dim_x = L;
	m_dim_y = L;
	m_Dim = m_dim_x * m_dim_y;
    NextNeigbor.resize(m_Dim);
}

Grid::Grid(size_t Lx, size_t Ly)
{
	m_dim_x = Lx;
	m_dim_y = Ly;
	m_Dim = m_dim_x * m_dim_y;
    NextNeigbor.resize(m_Dim);
}


/////////////////////////////GRIDFUNCTIONS/////////////////////////


/////////////////////////////VectorClass/////////////////////////
template<typename Val>
class Vector //1D Grid representation
{
    private:
    std::vector<Val> vec;

    Grid *grid; //get access to private of grid

    public:
    //Constructors
    Vector(Grid *grid);
    Vector(Grid *grid, std::mt19937 &gen);
		
    //Destructor
    ~Vector();
		
    //Operators
    Val& operator[](size_t index);
        
    //Functions
    size_t Dim() {return grid->Dim();}    
	
    void print();

};

//Constructors
template<typename Val>
Vector<Val>::Vector(Grid *grid) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);

     for(int point = 0; point< DIM; ++point)
    {
        vec[point]=-1; //cold start spins up
    }

}

template<typename Val>
Vector<Val>::Vector(Grid *grid,std::mt19937 &gen) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);
    std::uniform_int_distribution<int> unidist(0,100);
    
     for(int point = 0; point< DIM; ++point)
    {
        double RN = unidist(gen);
        if(RN <= 50)
        vec[point]=-1; //cold start spins up
        else
        vec[point]=1;
    }

}

//Destructor
template<typename Val>
Vector<Val>::~Vector()
{}

//operators
template<typename Val>
Val& Vector<Val>::operator[](size_t index)
{
	return vec[index];
}

/////////////////////////////VECTORFUNCTIONS/////////////////////////
template<typename Val>
void Vector<Val>::print() //prints whole vector
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
double ED(Vector<Val> &Configuration, Grid *g) //Energy density
{
    double Energy = 0;

    for(int point = 0; point< Configuration.Dim(); ++point)
    {
        Energy -= g->getB() * Configuration[point];

        for(int position = 0; position < 4; ++position)
        {
            Energy += -1*g->getJ() * Configuration[point] * Configuration[g->getNN(point,position)];
        }
    }
    return 0.5*Energy/g->Dim();
}

template<typename Val>
double ED_v(std::vector<Val> &Configuration, Grid *g) //Energy density
{
    double Energy = 0;

    for(int point = 0; point< Configuration.size(); ++point)
    {
        Energy -= g->getB() * Configuration[point];

        for(int position = 0; position < 4; ++position)
        {
            Energy += -1*g->getJ() * Configuration[point] * Configuration[g->getNN(point,position)];
        }
    }
    return 0.5*Energy/g->Dim();
}
//Magnetisationdensity
template<typename Val>
double MD(Vector<Val> &Configuration, Grid *g) //Energy density
{
    double Magnetization = 0;

    for(int point = 0; point< Configuration.Dim(); ++point)
    {
        Magnetization += Configuration[point];
    }
    return abs(Magnetization)/g->Dim();
}

template<typename Val>
double MD_v(std::vector<Val> &Configuration, Grid *g) //Energy density
{
    double Magnetization = 0;

    for(int point = 0; point< Configuration.size(); ++point)
    {
        Magnetization += Configuration[point];
    }
    return abs(Magnetization)/g->Dim();
}






/////////////////////////////LatticeClass/////////////////////////
template<typename Val>
class Lattice //2D representation of Grid
{
    private:
    std::vector<std::vector<Val>> Lat;

    Grid *grid;
    
    public:
    Lattice(Grid *grid);
    
    ~Lattice();

    Val& operator()(size_t i, size_t j);

    size_t DimX() {return grid->DimX();}
    size_t DimY() {return grid->DimY();}
    
    void print();
};

//Constructors
template<typename Val>
Lattice<Val>::Lattice(Grid *grid) :grid(grid)
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
Lattice<Val>::~Lattice()
{}

//Operators
template<typename Val>
Val& Lattice<Val>::operator()(size_t x, size_t y)
{
	return Lat[x][y];
}
/////////////////////////////LATTICEFUNCTIONS/////////////////////////
template<typename Val>
void Lattice<Val>::print()
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
Vector<Val> transform(Lattice<Val> &lat, Grid *grid)
{
	Vector<Val> vec(grid);

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
Lattice<Val> transform(Vector<Val> &vec, Grid *grid) //transforms vec to lat
{
	Lattice<Val> lat(grid);

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
void NN(Grid *g)
{
    //NextNeigbor = Vector(top,bottom,left,right)

    Lattice<Val> lat(g); //lat(x,y)

    int x,y = 0;

    const auto Lx = g->DimX();
    const auto Ly = g->DimY();

    for(int point = 0; point < g->Dim(); ++point)
    {
        x = point/Ly; //X position in Grid (up to down)
        y = point % Ly;    //Y position in Grid (left to right)

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


//Metropolis//////////////////////////////////////////////////

//Delta Energy
template<typename Val>
void q(double &q,Vector<Val> &Configuration, Grid *g, int &point)
{
    for(int i = 0; i<4; ++i)
    {q += Configuration[g->getNN(point,i)];}

    q *= Configuration[point];
}

//Spinflip
template<typename Val>
void Spinflip(Vector<Val> &Configuration,int point)
{
    Configuration[point] *= -1;
}


template<typename Val>
void Metropolis(Vector<Val> &Configuration, Grid *g, int Iterations, std::mt19937 &gen)
{
    int MarkovTime = 0;
/*
    const char* Datei = "IsingE100.txt";
    const char* Datei2 = "IsingM100.txt";
    FILE * handle = fopen(Datei, "w");
    FILE * handle2 = fopen(Datei2, "w");
*/
    int acceptance = 0;
    int proposals = 0;

    double exponential = 0;

    double q2 = exp(-4*g->getBeta()*g->getJ()); //"lookuptable"
    double q4 = exp(-8*g->getBeta()*g->getJ());

    double VarianceE = 0;
    double VarianceM = 0;

    vector<double> ED_vector;
    vector<double> MD_vector;

    vector<double> ED_vector_squared;
    vector<double> MD_vector_squared;


    double MDavg = 0;
    double EDavg = 0;
    
    
    //Markov Process
    for(int i=0; i< Iterations; ++i)
    {
        //fprintf(handle, "%d ",MarkovTime);
        //fprintf(handle2, "%d ",MarkovTime);

        //Getting each point of Configuration
        for(int point = 0; point < Configuration.Dim(); ++point)
        {
            double Q  = 0;
            q(Q,Configuration, g, point);
            Spinflip(Configuration,point); //new proposal

            proposals += 1;

            if(Q <= 0)
            {acceptance += 1;}

            if(Q > 0)
            {
                if(g->getB() == 0)
                {
                    if(Q == 2)
                    {exponential = q2;}
                    if(Q == 4)
                    {exponential = q4;}
                }

                else
                {
                    double Bexpo = exp(2*g->getB()*g->getBeta() *Configuration[point]);
                    if(Q == 2)
                    {exponential = q2 * Bexpo;}
                    if(Q == 4)
                    {exponential = q4 * Bexpo;}
                }

                double RandomNumber = 0;
                RNG_uni(RandomNumber,gen);
                //RandomNumber = RandomNumber/100;

                double Rho = exponential;

                if(RandomNumber > Rho)
                {Spinflip(Configuration,point);}

                else
                {acceptance += 1;}//else would be accept the new config
            }
        }

        //fprintf(handle, "%lf\n",ED(Configuration,g));
        //fprintf(handle2, "%lf\n",MD(Configuration,g));
        MarkovTime += 1;

        if(MarkovTime > 5000) //Expectationvalue of E and M after equilibration and after each Markov Chain
        {
            ED_vector.push_back(ED(Configuration,g));
            MD_vector.push_back(MD(Configuration,g));
            //std::cout << "lol" << std::endl;

            ED_vector_squared.push_back(pow(ED(Configuration,g),2));
            MD_vector_squared.push_back(pow(MD(Configuration,g),2));
        }
    }
    //std::cout << "lal" << std::endl;

    MDavg = Mean(MD_vector);
    EDavg = Mean(ED_vector);

    vector<double> C_eff;
    vector<double> X_eff;

    for(int m = 0; m < ED_vector.size(); ++m)
    {
        C_eff.push_back(pow(g->getBeta(),2)* Configuration.Dim() * (ED_vector_squared[m] - 2 * EDavg * ED_vector[m]));
        X_eff.push_back(g->getBeta() * Configuration.Dim() * (MD_vector_squared[m] - 2 * MDavg * MD_vector[m]));
    }    
  

    double tau_E = tau_int(ED_vector);
    double tau_E_squared = tau_int(ED_vector_squared);

    double tau_M = tau_int(MD_vector);
    double tau_M_squared = tau_int(MD_vector_squared);

    double tau_C = tau_int(C_eff);
    double tau_X = tau_int(X_eff);

    double sigma_E = auto_std_err_prim(ED_vector.size(),tau_E, ED_vector);
    double sigma_M = auto_std_err_prim(MD_vector.size(),tau_M, MD_vector);
    double sigma_C = auto_std_err_prim(C_eff.size(),tau_C, C_eff);
    double sigma_X = auto_std_err_prim(X_eff.size(),tau_X, X_eff);

    double C = pow(g->getBeta(),2)* Configuration.Dim() * svar(ED_vector);
    double X = g->getBeta() * Configuration.Dim() * svar(MD_vector);

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
    std::cout << g->getT() << " " <<tau_E << " " << tau_M << " " << tau_C <<  " " << tau_X << std::endl;


}

//WOLFF algorithm//////////////////////////////////////////////////

template<typename Val>
void Wolff(Vector<Val> &Configuration, Grid *g, int Iterations, std::mt19937 &gen)
{
    int sweep = 0;
    
    std::uniform_int_distribution<int> unidist(0,Configuration.Dim()-1);
    
    //accept probability
    double P_add = 1 - exp(-2*g->getBeta());
    
    vector<double> ED1_vector;
    vector<double> MD1_vector;
    vector<double> CS_vector;

    vector<double> ED1_vector_squared;
    vector<double> MD1_vector_squared;
    
    for(int i = 0; i < Iterations; ++i)
    {   
        std::vector<Val> Cluster;
        
        //choose random site
        int randomsite = unidist(gen);
        Cluster.push_back(randomsite);
        
        
        int new_entries = 1;
        int actual = 0;
        
        int spin = Configuration[randomsite];

        Spinflip(Configuration,randomsite);
        
        while(actual < new_entries)
        {
            randomsite = Cluster[actual];

            //go through the neigbors
            for(int neighbor = 0; neighbor < 4; ++neighbor)
            {
                //accept reject for all 4 neigbors
                //-->if spin is the same, accept reject step
                if(Configuration[g->getNN(randomsite,neighbor)] == spin)
                {
                    //Random number generator
                    double RN = 0;
                    RNG_uni(RN,gen);
                    //RN = RN/100;
           
                    //add to cluster if RN<=P_add
                    if(RN <= P_add)
                    {
                        int nn = g->getNN(randomsite,neighbor);
                        Spinflip(Configuration, nn);
                        
                        Cluster.push_back(g->getNN(randomsite,neighbor));
                        new_entries += 1;
                    }
                }
            }
            
            actual += 1;

        }
    sweep += 1;

    if(sweep > 5000) //Expectationvalue of E and M after equilibration and after each Markov Chain
        {
            ED1_vector.push_back(ED(Configuration,g));
            MD1_vector.push_back(MD(Configuration,g));

            ED1_vector_squared.push_back(pow(ED(Configuration,g),2));
            MD1_vector_squared.push_back(pow(MD(Configuration,g),2));

            CS_vector.push_back(Cluster.size()/(double) Configuration.Dim());
        }

    }
    double MDavg1 = Mean(MD1_vector);
    double EDavg1 = Mean(ED1_vector);

    vector<double> C_eff1;
    vector<double> X_eff1;

    for(int m = 0; m < ED1_vector.size(); ++m)
    {
        C_eff1.push_back(pow(g->getBeta(),2)* Configuration.Dim() * (ED1_vector_squared[m] - 2 * EDavg1 * ED1_vector[m]));
        X_eff1.push_back(g->getBeta() * Configuration.Dim() * (MD1_vector_squared[m] - 2 * MDavg1 * MD1_vector[m]));
    }

    double tau_E1 = tau_int(ED1_vector);
    double tau_E_squared1 = tau_int(ED1_vector_squared);

    double tau_M1 = tau_int(MD1_vector);
    double tau_M_squared1 = tau_int(MD1_vector_squared);

    double tau_C1 = tau_int(C_eff1);
    double tau_X1 = tau_int(X_eff1);
    double tau_CS = tau_int(CS_vector);

    double sigma_E1 = auto_std_err_prim(ED1_vector.size(),tau_E1, ED1_vector);
    double sigma_M1 = auto_std_err_prim(MD1_vector.size(),tau_M1, MD1_vector);
    double sigma_C1 = auto_std_err_prim(C_eff1.size(),tau_C1, C_eff1);
    double sigma_X1 = auto_std_err_prim(X_eff1.size(),tau_X1, X_eff1);

    double C1 = pow(g->getBeta(),2)* Configuration.Dim() * svar(ED1_vector);
    double X1 = g->getBeta() * Configuration.Dim() * svar(MD1_vector);
    double CS = Mean(CS_vector);
    double CS_error = auto_std_err_prim(CS_vector.size(),tau_CS, CS_vector);

    //std::cout << tau_C << " " << tau_X << std::endl;
    //std::cout << "tau C: "<<auto_std_err_prim(C_eff.size(),tau_C, C_eff) << std::endl;
    //std::cout << "tau X: "<<auto_std_err_prim(X_eff.size(),tau_X, X_eff) << std::endl;
    //std::cout << "tau MD: "<< auto_std_err_prim(MD_vector.size(),tau_M, MD_vector) << std::endl;
    //std::cout << "tau MD: "<< auto_std_err_prim(ED_vector.size(),tau_M, ED_vector) << std::endl;
    //std::cout << Blocking(MD_vector,20, g->getBeta() * Configuration.Dim()) << std::endl;
    //std::cout << Blocking(MD_vector,20, g->getBeta() * Configuration.Dim()) << std::endl;

    //std::cout << Bootstrap(MD_vector,tau_M, 1000, gen, g->getBeta()  * Configuration.Dim()) << std::endl;
    //std::cout << Bootstrap(MD_vector,tau_X, 1000, gen, g->getBeta()  * Configuration.Dim()) << std::endl;
    //std::cout << g->getT() << " " << EDavg1 << " " << MDavg1 << " " << C1 <<  " " << X1 << std::endl;
    //std::cout << g->getT() << " " <<sigma_E1 << " " << sigma_M1 << " " << sigma_C1 <<  " " << sigma_X1 << std::endl;
    //std::cout << g->getT() << " " <<tau_E1 << " " << tau_M1 << " " << tau_C1 <<  " " << tau_X1 << std::endl;
    std::cout << g->getT() <<" "<< CS << " " << CS_error << " " << tau_CS << std:: endl;

        
}
#endif //Grid_H
