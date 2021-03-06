#include"initialize.h"
#include <trng/yarn5s.hpp>
#include <trng/uniform01_dist.hpp>

void init_sim(System::simulation& sim)
{
    double velsum[sim.n_types][sim.n_dimensions];
    double vel2sum[sim.n_types];
    double scale[sim.n_types];

    #pragma omp target teams distribute parallel for
    for(int i = 0; i < sim.n_types;++i)
    {
        vel2sum[i]=0;
        scale[i]=0;
        for(int k = 0;k < sim.n_dimensions;++k)
        {
            velsum[i][k]=0;
        }
        #pragma omp parallel for
        for(int j = 0; j < sim.n_particles[i];++j)
        {
            #pragma omp parallel for
            for(int k = 0; k < sim.n_dimensions;++k)
            {
                trng::yarn5s R;

                int size = omp_get_num_threads();
                int rank = omp_get_thread_num();
                R.split(size,rank);

                trng::uniform01_dist<> unif;

                sim.position[i][j][k] = unif(R)*sim.box_size_limits[k];
                sim.velocity[i][j][k] = (unif(R) - 0.5);
                velsum[i][k]+=sim.velocity[i][j][k];
                vel2sum[i]+=pow(sim.velocity[i][j][k],2);
            }
        }
        #pragma omp parallel for
        for(int k = 0; k < sim.n_dimensions;++k)
        {
            velsum[i][k]=velsum[i][k]/sim.n_particles[i];
        }
        vel2sum[i]=vel2sum[i]/sim.n_particles[i];
        scale[i]=sqrt(sim.n_dimensions*BOLTZ_SI*sim.temperature[i]/(sim.mass[i]*vel2sum[i]));

        for(int j = 0; j < sim.n_particles[i];++j)
        {
            for(int k = 0; k < sim.n_dimensions;++k)
            {
                sim.velocity[i][j][k]=(sim.velocity[i][j][k]-velsum[i][k])*scale[i];
            }
        }

    }
}


