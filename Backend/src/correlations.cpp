/** @file */
#include "correlations.h"

void initialize_correlations(System::simulation& sim)
{
	sim.velocity_initial = sim.velocity;
}

void correlate(System::simulation& sim)
{
	int s = sim.step;

	//Taking local sum
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.correlation_velocity[s][i] = 0;

		#pragma omp parallel for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			double vdot = 0;
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				vdot += sim.velocity[i][j][k]*sim.velocity_initial[i][j][k];
			}
			sim.correlation_velocity[s][i] += vdot;
		}
	}

	// Taking global sum
	#pragma omp parallel for reduction(+:sim.correlation_velocity[s][sim.n_types])
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.correlation_velocity[s][sim.n_types] += sim.correlation_velocity[s][i];
	}

	// Dividing by the number of particles
	#pragma omp parallel for
	for (int i = 0; i < n_types; ++i)
	{
		sim.correlation_velocity[s][i] /= sim.n_particles[i];
	}

	//Taking global average
	sim.correlation_velocity[s][n_types] /= sim.numpartot;
}