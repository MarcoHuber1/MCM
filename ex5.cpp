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
    Grid_2D.setT(1.2);
    Grid_2D.setBeta(1/(Grid_2D.getT()));

    //NextNeigbor table
    NN<double>(&Grid_2D, Theta); //Generating table for next neigbors

    const char* Datei = "thermali.txt";
    const char* Datei2 = "sth.txt";

    //Thermalisation:
    int t_HMC = 100000; double t_LF = 0.3;
    double stepsize = 0.01;
    //double n = 2;
    HMC(&Grid_2D,Theta,gen,t_HMC,t_LF,stepsize, Datei);
    //Metropolis_XY(Theta, &Grid_2D, t_HMC, gen);


/*
    t_HMC = 1000;
    std::cout << "Doing..." <<std::endl;
    for(stepsize = 0.01; stepsize < 0.425; stepsize += 0.025)
    {

        for(n = 1; n<30; n+=2.5)
        {

            HMC(&Grid_2D,Theta,gen,t_HMC,t_LF,stepsize,Datei2,n);
            //std::cout << t_LF << std::endl;

        }
        std::cout << std::endl;
    }
    std::cout << "Done" <<std::endl;
*/

/*
    t_HMC = 3000;
    t_LF = 0.3;
    for(double temp = 0.1; temp <= 4.0; temp+=0)
    {

        Grid_2D.setT(temp);
        Grid_2D.setBeta(1/temp);

        //Spin_vector<double> Theta(&Grid_2D);

        HMC(&Grid_2D,Theta,gen,t_HMC,t_LF, stepsize);
        //Metropolis_XY(Theta, &Grid_2D, 100000, gen);

        temp += 0.1;
    }
*/
}
