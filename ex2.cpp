#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{

	Grid lattice(3);
	Vector<int> vec(&lattice);

	//Grid::Vector<double> Vec2(4);
//	Grid::Vector<double> Vec2(3,4);
//	std::cout << Vec[0];
	for(int i = 0; i<9; ++i)
		std::cout << vec[i] << std::endl;

	//std::cout << vec[0];
//	Vec + Vec2;
}

