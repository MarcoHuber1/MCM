#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

int main()
{
    //Set up random number generator
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = 1111151338;
    std::mt19937 gen;
    gen.seed(seed);

    //Set up Grid
    Grid Grid_2D(64); //argument is lattice dimension
    Lattice<int> Lat(&Grid_2D);



    //Configuration vector
    Vector<int> Configuration(&Grid_2D);
    Lat = transform(Configuration,&Grid_2D); //2D representation
    //Lat.print();

    /////////////////parameters for Grid//////////////
    double kb = 1.3806488 * pow(10,-23);
    Grid_2D.setJ(1);
    Grid_2D.setk(1);
    Grid_2D.setB(0);
    Grid_2D.setT(2.0);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<int>(&Grid_2D); //Generating table for next neigbors
    //Metropolis(Configuration,&Grid_2D,10000,gen);
    
/////////////////////////////////////////////////////////////////////////////////////    
    for(double temp = 1; temp <= 4.0; temp+=0)
    {

        if(temp>2.1 && temp<2.4)
        {
            
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            //unsigned seed = 1111121338;
            std::mt19937 gen;
            gen.seed(seed);
            
            Grid_2D.setT(temp);
            Grid_2D.setBeta(1/(Grid_2D.getT()));
            Vector<int> Configuration(&Grid_2D);
            
            //Metropolis(Configuration, &Grid_2D, 10000, gen);
            Wolff(Configuration, &Grid_2D, 10000, gen);
            temp += 0.005;
        }
        else
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            //unsigned seed = 1111121338;
            std::mt19937 gen;
            gen.seed(seed);
    
            Grid_2D.setT(temp);
            Grid_2D.setBeta(1/(Grid_2D.getT()));

            Metropolis(Configuration, &Grid_2D, 10000, gen);
            //Wolff(Configuration, &Grid_2D, 100000, gen);
            temp += 0.1;
        }

    }

/*
   for(double temp = 1; temp <= 4.0; temp+=0)
    {
            Vector<int> Configuration(&Grid_2D);
            Grid_2D.setT(temp);
            Grid_2D.setBeta(1/temp);

            Metropolis(Configuration, &Grid_2D, 20000, gen);
            //Wolff(Configuration, &Grid_2D, 10000, gen);
            temp += 0.02;

    }
*/
}
