/** @file */ 
#include "integrate.h"

/*******************************************************************************
 * This function integrates the equation of motion for periodic boundary conditions
 * This function integrates using the Leapfrog Algorithm 
 *
 * @param sim Simulation being integrated over
 ******************************************************************************/
void integrate_verdet_periodic(System::simulation& sim){
	double dt = sim.timestep;
	sim.energy_total = sim.energy_potential;
	//Starting first part of thermostat
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.position[i][j][k] += dt*sim.velocity[i][j][k] + 0.5*dt*dt*sim.acceleration[i][j][k];
				sim.position[i][j][k] -= sim.box_size_limits[k]*std::floor(sim.position[i][j][k]/sim.box_size_limits[k]);
				sim.velocity[i][j][k] += dt*sim.acceleration[i][j][k]/2;
			}
		}
	}
	//Done with first part of thermostat
	//Calling interaction
	sim.energy_potential = 0;
	interact(sim);
	//Done
	//Starting second part of thermostat
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.energy_kinetic[i] = 0;
		#pragma omp parallel for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.velocity[i][j][k] += dt*sim.acceleration[i][j][k]/2;
				#pragma omp atomic
				sim.energy_kinetic[i] += sim.velocity[i][j][k]*sim.velocity[i][j][k];
			}
		}
		sim.energy_kinetic[i] *= sim.mass[i];
		sim.temperature[i] = sim.energy_kinetic[i]/(sim.n_particles[i]*BOLTZ_SI*sim.n_dimensions);
		sim.energy_total += sim.energy_kinetic[i];
	}
	sim.time+= dt;
	sim.state++;
}

/*******************************************************************************
 * This function integrates the equation of motion for rigid box conditions
 * This function integrates using the Leapfrog Algorithm 
 *
 * @param sim Simulation being integrated over
 ******************************************************************************/
void integrate_verdet_box(System::simulation& sim){
	double dt = sim.timestep;
	sim.energy_total = sim.energy_potential;
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.energy_kinetic[i] = 0;
		#pragma omp target teams distribute parallel for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.position[i][j][k] += dt*sim.velocity[i][j][k] + 0.5*dt*dt*sim.acceleration[i][j][k];
				sim.velocity[i][j][k] += dt*sim.acceleration[i][j][k];
				if(sim.position[i][j][k]<0 || sim.position[i][j][k] > sim.box_size_limits[k]){
					int num_bounce = (int)(sim.position[i][j][k]/sim.box_size_limits[k]);
					double l = std::abs(sim.position[i][j][k] - sim.box_size_limits[k]*num_bounce);
					//Getting rid of the signs
					num_bounce *= ((num_bounce < 0)*-1 + (num_bounce>0));
					//Done
					//Implementing reflection
					sim.position[i][j][k] = (num_bounce%2 == 0)*l + (num_bounce%2 == 1)*(sim.box_size_limits[k]-l);
					sim.velocity[i][j][k] *= ((num_bounce%2==0) + (num_bounce%2 == 1)*(-1));
					//Done
				}
				#pragma omp atomic
				sim.energy_kinetic[i] += sim.velocity[i][j][k]*sim.velocity[i][j][k];
			}
		}
		sim.energy_kinetic[i] *= sim.mass[i];
		sim.temperature[i] = sim.energy_kinetic[i]/(sim.n_dimensions*sim.n_particles[i]*BOLTZ_SI);
		sim.energy_total += sim.energy_kinetic[i];
	}
	sim.time+=dt;
	sim.state++;
}