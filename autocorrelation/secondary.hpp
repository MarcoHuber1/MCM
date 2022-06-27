#ifndef autocorrelation_secondary_H
#define autocorrelation_secondary_H

//Error propagation//////////////////////////////////////////////////

//Blocking method//////////////////////////////////////////////////
template<typename Val>
double Blocking(std::vector<Val> Density, int Blocksize, double a)
{
    std::vector<Val> Q_block;
    
    int N = Density.size();
    
    //Generating blocks and making avg for each block, storing in Q_block
    for(int block = 0; block < Blocksize; ++ block)
    {
        //Sample
        std::vector<Val> Q;
        Q.resize(N/Blocksize);
    
        for(int i = 0; i< N/Blocksize; ++i)
        {
            Q[i] = Density[(block * N/Blocksize)+i];
        }
        Q_block.push_back(a * svar(Q));
    }
    
    //Compute block average for blockvariance
    double Blockavg = 0;
    for(int i = 0; i<Blocksize; ++i)
    {
        Blockavg += Q_block[i];
    }
    Blockavg = Blockavg / Blocksize;

    //Computing Block variance
    double BlockVariance = 0;
    for(int i = 0; i<Blocksize; ++i)
    {
        BlockVariance += pow(Q_block[i] - Blockavg,2);
    }
    BlockVariance = BlockVariance/Blocksize;

    //return std deviation
    return sqrt(BlockVariance)/sqrt(Blocksize);
}

//Bootstrap//////////////////////////////////////////////////
template<typename Val>
double Bootstrap(std::vector<Val> Density, double tau, int M, std::mt19937 gen, double a)
{
    int N = Density.size();
    int N_indep = N/(2*tau);
    std::uniform_int_distribution<int> unidist(0,N_indep);
    
    //Generate independent Sample with autocorr.time
    std::vector<Val> Sample;
    for(int i = 0; i< N_indep; ++i)
    {
        Sample.push_back(Density[i*2*tau]);
    }
    
    std::vector<Val> Q_boot;
    
    //generating M bootstrap samples Q
    for(int j = 0; j < M; ++j)
    {
        std::vector<Val> Q;
        Q.resize(N_indep);
        
        //giving each sample random numbers of density
        for(int i = 0; i < N_indep; ++i)
        {
            int randomplace = unidist(gen);
            Q[i] = Sample[randomplace];
        }
        Q_boot.push_back(svar(Q) * a); //give back mean of each bootstrap sample
    }
    
    //Compute Bootstrap average
    double Bootavg = 0;
    for(int m = 0; m<M; ++m)
    {
        Bootavg += Q_boot[m];
    }
    Bootavg = Bootavg/M;
    
    //Computing Boot std
    double BootDeviation = 0;
    for(int n = 0; n<M; ++n)
    {
        BootDeviation += pow(Q_boot[n] - Bootavg,2);
    }
    return sqrt(BootDeviation/M);
}
#endif // autocorrelation_secondary_H
