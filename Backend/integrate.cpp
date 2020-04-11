#include "integrate.h"

void integrate_verdet_periodic(simulation* sim){
	dt = sim.timestep;
	sim.energy_total = sim.energy_potential;
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.energy_kinetic[i] = 0;
		#pragma omp target teams distribute for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < n_dimensions; ++k)
			{
				sim.position[i][j][k] += dt*sim.velocity[i][j][k] + 0.5*dt*dt*sim.acceleration[i][j][k];
				sim.position[i][j][k] -= sim.box_size_limits[k]*std::floor(sim.position[i][j][k]/sim.box_size_limits[k]);
				sim.velocity[i][j][k] += dt*sim.acceleration[i][j][k];
				#pragma omp atomic
				{
					sim.energy_kinetic[i] += sim.velocity[i][j][k]*sim.velocity[i][j][k];
				}
			}
		}
		sim.energy_kinetic[i] *= sim.mass[i];
		sim.temperature[i] = sim.energy_kinetic[i]/(sim.n_dimensions*sim.n_particles*BOLTZ_SI);
		sim.energy_total += sim.energy_kinetic[i];
	}
}