#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "interaction.h"

double distance(System::simulation& sim, int type1, int n1, int type2, int n2)
{
	double dist = 0;
	if(sim.periodic_boundary){
		for (int i = 0; i < sim.n_dimensions; ++i)
		{
			double temp = abs(sim.position[type1][n1][i] - sim.position[type2][n2][i]);
			double temp2 = std::min(sim.box_size_limits[i] - temp,temp);
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
	return dist;
}

void interact(System::simulation& sim){
	sim.energy_potential = 0;
	for (int i = 0; i < sim.n_types; ++i)
	{
		for (int j = 0; j < sim.n_types; ++j)
			{
				sim.interaction[i][j](sim,i,j);
			}		
	}
}

void free_particles(System::simulation& sim, int type1, int type2)
{
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type1]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			sim.acceleration[type1][j][k] = 0;
		}
	}
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type2]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			sim.acceleration[type2][j][k] = 0;
		}
	}
}

void lj_periodic(System::simulation& sim){
	
}