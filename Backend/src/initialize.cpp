#include"initialize.h"

void init_sim(System::simulation& sim)
{
    double velsum[sim.n_types][sim.n_dimensions];
    double vel2sum[sim.n_types];
    double scale[sim.n_types];

    //Initializing the randomization
    std:: random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> unif(0, 1);
    //Randomization Initialized


    for(int i = 0; i < sim.n_types;++i)
    {
        vel2sum[i]=0;
        scale[i]=0;
        for(int k = 0;k < sim.n_dimensions;++k)
        {
            velsum[i][k]=0;
        }
        for(int j = 0; j < sim.n_particles[i];++j)
        {
            for(int k = 0; k < sim.n_dimensions;++k)
            {
                sim.position[i][j][k] = unif(gen)*sim.box_size_limits[k];
                sim.velocity[i][j][k] = (unif(gen) - 0.5);
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


