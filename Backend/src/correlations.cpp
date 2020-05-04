/** @file */
#include "correlations.h"

void initialize_correlations(System::simulation& sim)
{
	sim.velocity_initial = sim.velocity_com;
}

void correlate(System::simulation& sim)
{
	int s = sim.state;

	//Taking local sum
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.correlation_velocity[s][i] = 0;

		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			double vdot = 0;
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				vdot += sim.velocity_com[i][j][k]*sim.velocity_initial[i][j][k];
			}
			sim.correlation_velocity[s][i] += vdot;
		}
	}

	// Taking global sum
	double temp = sim.correlation_velocity[s][sim.n_types];
	#pragma omp parallel for reduction(+:temp)
	for (int i = 0; i < sim.n_types; ++i)
	{
		temp += sim.correlation_velocity[s][i];
	}
	sim.correlation_velocity[s][sim.n_types] = temp;

	// Dividing by the number of particles
	#pragma omp parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.correlation_velocity[s][i] /= sim.n_molecules[i];
	}

	//Taking global average
	sim.correlation_velocity[s][sim.n_types] /= sim.numpartot;
}