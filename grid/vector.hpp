#ifndef Grid_H
#define Grid_H

template<>
class Grid
{

	private:
     	//2D
	const int m_dim_x = 10;
	const int m_dim_y = 10;

	//1D Linearisation
	std::vector<double> Point(m_dim_x * m_dim_y,0);

	public:

	//Constructors
	Grid() = default;
	Grid(int L): m_dim_x(L), m_dim_y(L) {}
	Grid(int Lx, int Ly): m_dim_x(Lx), m_dim_y(Ly) {}

	//Get dimensions
	static constexpr auto dimX() const
	{
		return m_dim_x;
	}
	
	static constexpr auto dimY() const
	{
		return m_dim_y;
	}

};

#endif //Grid_H

