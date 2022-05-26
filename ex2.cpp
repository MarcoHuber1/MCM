#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{

	Grid lattice(3);
	Vector<int> vec(&lattice);
	
	Grid lattice2(3,4);
	Lattice<int> lat(&lattice2);
	lat(0,0) = 1;
	lat(0,1) = 2;
	lat(1,0) = 3;
	lat.print();

	Vector<int> trans(&lattice2); 
	trans = transform(lat,&lattice2);
	trans.print();

	Lattice<int> trans2(&lattice2);
	trans2 = transform(trans,&lattice2);
	trans2.print();

}

