#include"../include/client.h"

void one_step_md(simulation* sim)
{
    interact(sim);
    integrate_verdet_periodic(sim);
    write_traj(sim);
}