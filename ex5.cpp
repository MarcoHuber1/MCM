#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>
#include<ctime>

int main()
{
    //Set up random number generator
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = 1111121338;
    std::mt19937 gen;
    gen.seed(seed);

    //Set up Grid
    Grid_XY Grid_2D(16); //argument is lattice dimension
    Lattice_XY<double> Lat(&Grid_2D);

    //Configuration vector
    Spin_vector<double> Theta(&Grid_2D);

    Lat = transform(Theta,&Grid_2D); //2D representation
    //Lat.print();
    /////////////////parameters//////////////
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(0.3);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<double>(&Grid_2D, Theta); //Generating table for next neigbors
    double n = 1;
    int t_HMC = 100; double t_LF = 100;
    //HMC(&Grid_2D,Theta,gen,t_HMC,t_LF,n);
    //Metropolis_XY(Theta, &Grid_2D, t_HMC, gen);
/*
    for(double n=0; n<50; n+=0.1)
    {
        HMC(&Grid_2D,Theta,gen,t_HMC,t_LF,n);
        t_LF +=n;

    }
*/

    std::cout << "dt: 0.22" << " n: 15 "   << std::endl;
    for(double temp = 0.1; temp <= 4.0; temp+=0)
    {

        Grid_2D.setT(temp);
        Grid_2D.setBeta(1/temp);

        //Spin_vector<double> Theta(&Grid_2D);

        HMC(&Grid_2D,Theta,gen,t_HMC,t_LF, n);
        //Metropolis_XY(Theta, &Grid_2D, 10000, gen);
        temp += 0.1;
    }

}
