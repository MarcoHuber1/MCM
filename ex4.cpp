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
    Grid Grid_2D(64); //argument is lattice dimension
    Lattice<int> Lat(&Grid_2D);



    //Configuration vector
    Vector<int> Configuration(&Grid_2D, gen);
    Lat = transform(Configuration,&Grid_2D); //2D representation
    //Lat.print();

    /////////////////parameters//////////////
    double kb = 1.3806488 * pow(10,-23);
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(1.5);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<int>(&Grid_2D, Configuration); //Generating table for next neigbors
    //energydensity and magnetization

    double energy = ED<int>(Configuration, &Grid_2D);
    //std::cout << energy << std::endl;

    double magnetization = MD<int>(Configuration, &Grid_2D);


    //std::cout << magnetization << std::endl;

    for(double temp = 1; temp <= 4.0; temp+=0)
    {

        if(temp>2.1 && temp<2.4)
        {
            Grid_2D.setT(temp);
            Grid_2D.setBeta(1/(Grid_2D.getT()));

            Metropolis(Configuration, &Grid_2D, 100000, gen);
            //Wolff(Configuration, &Grid_2D, 10000, gen);
            temp += 0.02;
        }
        else
        {
            Grid_2D.setT(temp);
            Grid_2D.setBeta(1/(Grid_2D.getT()));

            Metropolis(Configuration, &Grid_2D, 100000, gen);
            //Wolff(Configuration, &Grid_2D, 10000, gen);
            temp += 0.1;
        }

    }


}
