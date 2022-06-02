#include<iostream>
#include<statistics.hpp>
#include<grid.hpp>

double Energy(std::vector<std::array<int,4>> NN, Vector<int> Gridvector,double J,double B);
double Z(Vector<int> Gridvector, double beta, double energy);

int main()
{
    /////////////////parameters//////////////
    double J = 1;
    double k = 1;
    double B = 1;
    double beta = -1;

    ////////////////Grid and Gridvector//////
    using table = std::vector<std::array<int,4>>;

	Grid Grid_2D(10);
	Vector<int> Gridvector(&Grid_2D);

    NN<int>(&Grid_2D); //Generating table for next neigbors
    double energy = ED<int>(&Gridvector, &Grid_2D);
    std::cout << energy << std::endl;
    
    double Zp = Z(Gridvector, beta, energy);
    std::cout << Zp << std::endl;
    
    
    
	

}


double Z(Vector<int> Gridvector, double beta, double energy)
{
    double Z = 0;
    for(int point = 0; point< Gridvector.Dim(); ++point)
    {
        Z += exp(-beta * energy);
    }
    return Z;
}
        
    
