#ifndef Grid_H
#define Grid_H

template<>
class Grid
{

	private:
     	//2D
	const int m_dim_x = 10;
	const int m_dim_y = 10;
	const int Dim = m_dim_x *m_dim_y;
	double PointXY[m_dim_x][m_dim_y];

	//1D Linearisation
	std::vector<double> Point(Dim,0);

	public:

	//Constructors
	Grid() = default;
	Grid(int L): m_dim_x(L), m_dim_y(L) {}
	Grid(int Lx, int Ly): m_dim_x(Lx), m_dim_y(Ly) {}

	//Gridmatrix:
	static constexpr auto GetSpin(int x, int y) const
	{
		return Position[x][y];
	}
	auto ChangeSpin(int x, int y, double value)
	{
		Position[x][y] = value;
	}
	
	//Reshaping
	//M->V
	auto reshapeMV const
	{
		for(
	//Get dimensions of grid
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

