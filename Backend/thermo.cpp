#include"thermo.h"

void trans_ke(simulation& sim)
{
    double vel2sum[sim.ntypes];
    for(int i = 0; i < sim.n_types;++i)
    {
        vel2sum[i]=0;
        for(int j = 0; j < sim.n_particles[i];++j)
        {
            for(int k = 0; k < sim.n_dimensions;++k)
            {
                vel2sum[i]+=pow(sim.velocity[i][j][k],2);
            }
        }
        vel2sum[i] *= 0.5*sim.mass[i];
        std::cout<<"KE of type "<<i<<" is "<<vel2sum[i]<<std::endl;
    
    }
}