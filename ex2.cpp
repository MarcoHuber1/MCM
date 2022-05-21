#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{
	Grid::Vector<double> Vec(4);
	//Grid::Vector<double> Vec2(3,4);
	std::cout << Vec[0];
	Vec[0] = 7;
	std::cout << Vec[0];
}

