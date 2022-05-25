#ifndef Grid_H
#define Grid_H
#include<statistics.hpp>
#include<cstddef>

class Grid
{
	private:
	size_t m_dim_x = 10;
	size_t m_dim_y = 10;

	public:	
	template<typename Val>	
	class Vector
	{
		private:
		Val *vec = nullptr;

		Grid *grid; //get access to private of grid
		size_t m_Dim = grid->m_dim_x * grid->m_dim_y; //Size of Vector
		
				
		public:
		//Constructors
		Vector(Grid *grid): grid(*grid) {}; 

		Vector() = default;	
		Vector(size_t L); //symm. grid
        	Vector(size_t Lx, size_t Ly); //not symm. gird
        	Vector(const Vector &orig);
		
		//Destructor
		~Vector();
		
		//Operators
		Val& operator[](size_t index) const;

		//Functions
		
	template<typename Val>	
	class Lattice
	{
		private:
		Val *mX = nullptr;
		Val *mY = nullptr;
		Val **matrix = nullptr;
		Grid *grid;
		Vector *vec;

		public:
		//Constructors
		Lattice(Vector *vec): vec(*vec) {};
		Lattice(Grid *grid): grid(*grid) {};
		
		Lattice() = default;
		Lattice(size_t L);
		Lattice(size_t Lx, size_t Ly);
		Lattice(const Lattice &orig);

	};

	};
};

/////////////////////Constructors//////////////////////////

template<typename Val>
Grid::Vector<Val>::Vector(size_t L)
{
	grid->m_dim_x = L;
	grid->m_dim_y = L;
	m_Dim = grid->m_dim_x * grid->m_dim_y;

	vec = new Val[m_Dim]{};
}

template<typename Val>
Grid::Vector<Val>::Vector(size_t Lx, size_t Ly)
{
        grid->m_dim_x = Lx;
        grid->m_dim_y = Ly;
        m_Dim = grid->m_dim_x * grid->m_dim_y;

        vec = new Val[m_Dim]{};
}

//Copy constr.
template<typename Val>
Grid::Vector<Val>::Vector(const Vector &orig): Vector(orig.m_Dim)
{
	delete[] vec;
	m_Dim = orig.m_Dim;
	vec = new Val[m_Dim]{};
	std::copy(orig.vec, orig.vec + m_Dim, vec);
}

template<typename Val>
Grd::Lattice<Val>::Lattice(size_t L)
{
	grid->m_dim_x = L;
	grid->m_dim_y = L;

	mX = new Val[m_dim_x];
	mY = new Val[m_dim_y];

	matrix = std::vector< 


}
/////////////////////Destructor///////////////////////////
template<typename Val>
Grid::Vector<Val>::~Vector()
{delete[] vec;}

////////////////////Operators////////////////////////////
template<typename Val>
Val& Grid::Vector<Val>::operator[](size_t index) const
{
	if(index>m_Dim)
	{std::cout << "Index out of range" << std::endl;}
	
	return vec[index];
}


#endif //Grid_H

