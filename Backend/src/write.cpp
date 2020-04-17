#include"write.h"

// for now this just dumps positions and velocities on the screen on the screen. later we can use it to write onto a file.
void write_traj(simulation& sim)
{
    for(int i = 0; i < sim.n_types;++i)
    {
        for(int j = 0; j < sim.n_particles[i];++j)
        {
            std::cout<<sim.state<<","<<i<<","<<j<<",";
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.position[i][j][k]<<",";
            }
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.velocity[i][j][k]<<",";
            }
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.acceleration[i][j][k]<<",";
            }
            std::cout<<std::endl;
        }   
    
    }
}