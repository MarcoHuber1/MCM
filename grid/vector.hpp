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
		size_t m_Dim = m_dim_x *m_dim_y;
	

		public:
		Vector() = default;
		Vector(size_t L); //symm. grid
        	Vector(size_t Lx, size_t Ly); //not symm. gird
        	Vector(const Vector &orig);

		~Vector();
	};
};

/////////////////////Constructors//////////////////////////

template<typename Val>
Grid::Vector<Val>::Vector(size_t L)
{
	m_dim_x = L;
	m_dim_y = L;
	m_Dim = m_dim_x * m_dim_y;

	vec = new Val[m_Dim]{};
}

template<typename Val>
Grid::Vector<Val>::Vector(size_t Lx, size_t Ly)
{
        m_dim_x = Lx;
        m_dim_y = Ly;
        m_Dim = m_dim_x * m_dim_y;

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

/////////////////////Destructor///////////////////////////



#endif //Grid_H

