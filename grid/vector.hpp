#ifndef Grid_H
#define Grid_H
#include<statistics.hpp>
#include<cstddef>


class Grid
{
	private:
	size_t m_dim_x = 10;
	size_t m_dim_y = 10;
	size_t m_Dim = m_dim_x * m_dim_y;

	public:	
	Grid(size_t Lx, size_t Ly);
	Grid(size_t L);
    
    //Functions
    size_t Dim() {return m_Dim;};
    size_t DimX() {return m_dim_x;};
    size_t DimY() {return m_dim_y;};
};

template<typename Val>
class Vector
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

template<typename Val>
class Lattice
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

/////////////////////Constructors//////////////////////////

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

template<typename Val>
Vector<Val>::Vector(Grid *grid) :grid(grid)
{
    auto DIM = grid->Dim();
    vec.resize(DIM);

    for(int i = 0; i<DIM; ++i)
    {vec[i] = 0;}
}

template<typename Val>
Lattice<Val>::Lattice(Grid *grid) :grid(grid)
{
    auto X = grid->DimX();
    auto Y = grid->DimY();
    
    Lat.resize(Y);
    for(int i=0; i<Y; ++i)
    {
        Lat[i].resize(X);
    }
 
}

/////////////////////Destructor///////////////////////////
template<typename Val>
Vector<Val>::~Vector()
{}

template<typename Val>
Lattice<Val>::~Lattice()
{}

////////////////////Operators////////////////////////////
template<typename Val>
Val& Vector<Val>::operator[](size_t index)
{
	return vec[index];
}

template<typename Val>
Val& Lattice<Val>::operator()(size_t x, size_t y)
{
	return Lat[x][y];
}

//////////////////Functions/////////////////////////////
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
}

template<typename Val>
void Vector<Val>::print()
{
	auto D = Dim();

	for(int d = 0; d<D; ++d)
	{
		std::cout << vec[d] << " ";
	}
	std::cout << "\n";
}

template<typename Val>
Lattice<Val> transform(Vector<Val> &vec, Grid *grid)
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
			vec[i*X+j] = lat(i,j);
		}
	}

	return vec;
}

#endif //Grid_H
