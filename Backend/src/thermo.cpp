#include"thermo.h"

void trans_ke(System::simulation& sim)
{
    double vel2sum[sim.n_types];
    #pragma omp target teams distribute parallel for
    for(int i = 0; i < sim.n_types;++i)
    {
        vel2sum[i]=0;

        #pragma omp parallel for
        for(int j = 0; j < sim.n_molecules[i];++j)
        {
            #pragma omp parallel for
            for(int k = 0; k < sim.n_dimensions;++k)
            {
                vel2sum[i]+=pow(sim.mol_state[i][j].velocity_com[k],2);
            }
        }
        vel2sum[i] *= 0.5*sim.mass[i];    
    }
    for (int i = 0; i < sim.n_types; ++i)
    {
        std::cout<<"KE of type "<<i<<" is "<<vel2sum[i]<<std::endl;
    }
}