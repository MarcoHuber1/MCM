#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{

	Grid lattice(3);
	Vector<int> vec(&lattice);
    
	Grid lattice2(3,8);
	Lattice<int> lat(&lattice2);
    lat.print();
	

	Vector<int> trans(&lattice2); 
	trans = transform(lat,&lattice2);
	trans.print();
     
	Lattice<int> trans2(&lattice2);
    
	trans2 = transform(trans,&lattice2);
    
	trans2.print();
   
    auto Next = NN<int>(&lattice2);
    std::cout << "top"<< Next[23][0] << std::endl;
    std::cout << "bot"<<Next[23][1] << std::endl;
    std::cout << "left"<<Next[23][2] << std::endl;
    std::cout << "right"<<Next[23][3] << std::endl;

}

