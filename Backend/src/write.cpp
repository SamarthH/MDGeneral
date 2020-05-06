#include"write.h"

// for now this just dumps positions and velocities on the screen on the screen. later we can use it to write onto a file.
void write_traj(System::simulation& sim)
{
    for(int i = 0; i < sim.n_types;++i)
    {
        for(int j = 0; j < sim.n_molecules[i];++j)
        {
            std::cout<<sim.state<<","<<i<<","<<j<<",";
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.mol_state[i][j].position_com[k]<<",";
            }
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.mol_state[i][j].velocity_com[k]<<",";
            }
            for(int k = 0; k < sim.n_dimensions;++k)
            {
               std::cout<<sim.mol_state[i][j].acceleration_com[k]<<",";
            }
            std::cout<<std::endl;
        }   
    
    }
}