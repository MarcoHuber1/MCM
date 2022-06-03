#ifndef Grid_H
#define Grid_H
#include<statistics.hpp>
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
    double mB = 1;
    double mbeta = -1;
    double mT = 0;

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
        vec[point]=1;
    }
    for(int point = 0; point< DIM/3; ++point)
    {
        vec[3*point]=-1;
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
        for(int position = 0; position < 4; ++position)
        {
            Energy += -g->getJ() * Configuration[point] * Configuration[g->getNN(point,position)] - g->getB() * Configuration[point];

            //std::cout << Configuration[g->getNN(point,position)] <<std::endl;
        }
    }
    return Energy/g->Dim();
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
void NN(Grid *g, Vector<Val> &Configuration)
{
    //NextNeigbor = Vector(top,bottom,left,right)

    Lattice<Val> lat(g); //lat(x,y)

    const auto Lx = g->DimX();
    const auto Ly = g->DimY();

    for(int point = 0; point < g->Dim(); ++point)
    {
            int x = (int)point/Ly; //X position in Grid (up to down)
            int y = point % Ly;    //Y position in Grid (left to right)

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


//Markov process//////////////////////////////////////////////////

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
void Spinflip(Vector<Val> &Configuration,int &point)
{
    Configuration[point] *= -1;
}


double lookuptable(double &Q,int &Spin, Grid *g)
{
    double exponential = 0;
    double Bexpo = 2*g->getB()*g->getBeta() *Spin;

    double q2 = exp(-4*g->getBeta()*g->getJ());
    double q4 = exp(-8*g->getBeta()*g->getJ());


    if(g->getB() == 0)
    {
        if(Q == 2)
        {exponential = q2;}
        if(Q == 4)
        {exponential = q4;}
    }
    else
    {
        if(Q == 2)
        {exponential = q2 *Bexpo;}
        if(Q == 4)
        {exponential = q4 *Bexpo;}
    }
    return exponential;

}

template<typename Val>
void Markov(Vector<Val> &Configuration, Grid *g, int Iterations, std::mt19937 &gen)
{
    int MarkovTime = 0;
    Vector<Val> rho(g);
    Lattice<Val> lat(g);
    lat = transform(Configuration, g);
    //lat.print();

    int acceptance = 0;
    int rejection = 0;
    int proposals = 0;
    //Markov iteration
    for(int i=0; i< Iterations; ++i)
    {
        //std::cout << "MT: "<< MarkovTime << std::endl;

        //Getting each point of Configuration
        for(int point = 0; point < Configuration.Dim(); ++point)
        {
            //lat = transform(Configuration, g);
            //lat.print();
            Spinflip(Configuration,point); //new proposal
            proposals +=1;
            double Q = 0;
            q(Q,Configuration, g, point);
            //std::cout << Q << std::endl;
            if(Q < 0)
                acceptance += 1;

            if(Q > 0)
            {
                double RandomNumber = 0;
                RNG_uni(RandomNumber,gen);
                RandomNumber = RandomNumber/100;
                double Rho = lookuptable(Q,Configuration[point],g);
                //std::cout << "Rho: "<< Rho <<std::endl;
                //std::cout << "RN: "<< RandomNumber <<std::endl;

                if(RandomNumber <= Rho)
                {
                    Spinflip(Configuration,point);
                    //std::cout << "rejected" << std::endl;
                    acceptance += 1;
                }
                else
                {
                    //std::cout << "accepted" << std::endl;
                    rejection += 1;
                }

            }
            //else would be accept the new config

        }
        MarkovTime += 1;
    }
    std::cout << (double)acceptance/(proposals) << std::endl;
}
#endif //Grid_H
