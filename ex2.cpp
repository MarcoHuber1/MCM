#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>


int main()
{
    //Set up Grid
    Grid Grid_2D(10); //argument is lattice dimension
    Lattice<int> Lat(&Grid_2D);


    //Configuration vector
	Vector<int> Configuration(&Grid_2D);
    Lat = transform(Configuration,&Grid_2D);
    Lat.print();

    /////////////////parameters//////////////
    double kb = 1.3806488 * pow(10,-23);
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(1*kb);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    using table = std::vector<std::array<int,4>>;
    NN<int>(&Grid_2D); //Generating table for next neigbors

    //energydensity and magnetization

    double energy = ED<int>(Configuration, &Grid_2D);
    std::cout << energy << std::endl;

    double magnetization = MD<int>(Configuration, &Grid_2D);
    std::cout << magnetization << std::endl;
    


    
    
    
	

}

        
    
