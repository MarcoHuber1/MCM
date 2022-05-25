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

/*		
	class Lattice
	{
		private:
		Val *mX = nullptr;
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
*/
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
	vec.resize(grid->GetGrid_Dim());
}

/////////////////////Destructor///////////////////////////
template<typename Val>
Vector<Val>::~Vector()
{}

////////////////////Operators////////////////////////////
template<typename Val>
Val& Vector<Val>::operator[](size_t index)
{
	return vec[index];
}

#endif //Grid_H
