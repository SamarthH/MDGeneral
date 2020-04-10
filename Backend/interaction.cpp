#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "interaction.h"

double distance(simulation* sim, int type1, int n1, int type2, int n2)
{
	double dist = 0;
	if(periodic_boundary){
		for (int i = 0; i < sim.n_dimensions; ++i)
		{
			double temp = abs(sim.position[type1][n1][i] - sim.position[type2][n2][i]);
			double temp2 = min(sim.box_size_limits[i] - temp,temp);
			dist+= temp2*temp2;
		}
		dist = std::sqrt(dist);
	}
	else{
		for (int i = 0; i < sim.n_dimensions; ++i)
		{
			double temp2 = abs(sim.position[type1][n1][i] - sim.position[type2][n2][i]);
			dist+= temp2*temp2;
		}
		dist = std::sqrt(dist);
	}
}

void free_particles_periodic(simulation* sim)
{
	double dt = sim.timestep;
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp target teams distribute for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < n_dimensions; ++k)
			{
				sim.position[i][j][k] = std::fmod(sim.position[i][j][k] + sim.velocity[i][j][k]*dt, box_size_limits[k]);
			}
		}
	}
}

void free_particles_box(simulation* sim)
{
	double dt = sim.timestep;
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp target teams distribute for
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < n_dimensions; ++k)
			{
				sim.position[i][j][k] = sim.position[i][j][k] + sim.velocity[i][j][k]*dt;
				while(sim.position[i][j][k] >= box_size_limits){
					sim
				}
			}
		}
	}
}