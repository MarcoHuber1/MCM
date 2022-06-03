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
    void setNN(int i, int j, double value){NextNeigbor[i][j] = value;};
    void resizeNN(int size){NextNeigbor.resize(m_Dim);}
    double getNN(int i, int j){return NextNeigbor[i][j];};




};

Grid::Grid(size_t L)
{
	m_dim_x = L;
	m_dim_y = L;
	m_Dim = m_dim_x * m_dim_y;
}

Grid::Grid(size_t Lx, size_t Ly)
{
	m_dim_x = Lx;
	m_dim_y = Ly;
	m_Dim = m_dim_x * m_dim_y;
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

    double Energy();
};

//Constructors
template<typename Val>
Vector<Val>::Vector(Grid *grid) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);

    for(int i = 0; i<DIM; ++i)
    {vec[i] = i;}
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
			Lat[i][j] = 0;
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
    g->resizeNN(g->Dim());
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
#endif //Grid_H
