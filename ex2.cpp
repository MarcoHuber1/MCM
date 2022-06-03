#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>


int main()
{
    //Set up random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen;
    gen.seed(seed);

    //Set up Grid
    Grid Grid_2D(4); //argument is lattice dimension
    Lattice<int> Lat(&Grid_2D);


    //Configuration vector
	Vector<int> Configuration(&Grid_2D);
    Lat = transform(Configuration,&Grid_2D); //2D representation
    //Lat.print();

    /////////////////parameters//////////////
    double kb = 1.3806488 * pow(10,-23);
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(2.6);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<int>(&Grid_2D, Configuration); //Generating table for next neigbors
    //energydensity and magnetization

    double energy = ED<int>(Configuration, &Grid_2D);
    std::cout << energy << std::endl;

    double magnetization = MD<int>(Configuration, &Grid_2D);
    std::cout << magnetization << std::endl;

    Markov(Configuration, &Grid_2D, 1, gen);



    
    
    
	

}

        
    
