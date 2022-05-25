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
    size_t GetGrid_Dim() {return m_Dim;};
    size_t GetGrid_Dimx() {return m_dim_x;};
    size_t GetGrid_Dimy() {return m_dim_y;};
};
template<typename Val>
class Vector
{
    private:
    std::vector<Val> vec;

    Grid *grid; //get access to private of grid
				
    public:
    //Constructors
    Vector<Val>(Grid *grid);
		
    //Destructor
    ~Vector();
		
    //Operators
    Val& operator[](size_t index);
        
    //Functions
};

template<typename Val>
class Lattice
{
    private:
    std::vector<std::vector<Val>> Lat;

    Grid *grid;
    
    public:
    Lattice<Val>(Grid *grid);
    
    ~Lattice();

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
    auto DIM = grid->GetGrid_Dim();
	vec.resize(DIM);
    for(int i = 0; i<DIM; ++i)
    {vec[i] = 1;}
}

template<typename Val>
Lattice<Val>::Lattice(Grid *grid) :grid(grid)
{
    auto dim_x = grid->GetGrid_Dimx();
    auto dim_y = grid->GetGrid_Dimy();
    
    Lat.resize(dim_y);
    for(int i=0; i<dim_y; ++i)
    {
        Lat[i].resize(dim_x);
    }
    

    for(int i = 0; i<dim_x; ++i)
    {
        Lat[i][0] = 1;
        Lat[0][i] = -1;
    
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

#endif //Grid_H
