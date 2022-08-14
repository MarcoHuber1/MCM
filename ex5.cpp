#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{
    //Set up random number generator
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = 1111121338;
    std::mt19937 gen;
    gen.seed(seed);

    //Set up Grid
    Grid_XY Grid_2D(5); //argument is lattice dimension
    Lattice_XY<double> Lat(&Grid_2D);

    //Configuration vector
    Spin_vector<double> Theta(&Grid_2D, gen);
    Spin_vector<double> Spin_x(&Grid_2D, gen);
    Spin_vector<double> Spin_y(&Grid_2D, gen);
    
    Lat = transform(Theta,&Grid_2D); //2D representation
    Lat.print();
    /////////////////parameters//////////////
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(2.3);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<double>(&Grid_2D, Theta); //Generating table for next neigbors
    std::cout << ED_XY(Theta,&Grid_2D) << std::endl;
    std::cout << MD_XY(Theta,&Grid_2D) << std::endl;

    std::cout << Theta[0] << std::endl;

    int t_HMC = 1; int t_LF = 2;
    HMC(&Grid_2D,Theta,gen,t_HMC,t_LF);
    std::cout << Theta[0] << std::endl;
    std::cout << ED_XY(Theta,&Grid_2D) << std::endl;
    std::cout << MD_XY(Theta,&Grid_2D) << std::endl;
    Lat = transform(Theta,&Grid_2D); //2D representation
    Lat.print();
    //energydensity and magnetization

    //double energy = ED<int>(Configuration, &Grid_2D);
    //std::cout << energy << std::endl;

    //double magnetization = MD<int>(Configuration, &Grid_2D);

    //Wolff(Configuration, &Grid_2D, 10000, gen);
    //std::cout << magnetization << std::endl;
/*
    for(double temp = 0.1; temp < 5.0; temp+=0.1)
    {
        Grid_2D.setT(temp);
        Grid_2D.setBeta(1/(Grid_2D.getT()));
        Markov(Configuration, &Grid_2D, 100000, gen);
    }
    */

}
