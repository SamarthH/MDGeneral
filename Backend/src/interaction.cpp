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
		for (int i = 0; i < (*sim).n_dimensions; ++i)
		{
			double temp = abs((*sim).position[type1][n1][i] - (*sim).position[type2][n2][i]);
			double temp2 = min((*sim).box_size_limits[i] - temp,temp);
			dist+= temp2*temp2;
		}
		dist = std::sqrt(dist);
	}
	else{
		for (int i = 0; i < (*sim).n_dimensions; ++i)
		{
			double temp2 = abs((*sim).position[type1][n1][i] - (*sim).position[type2][n2][i]);
			dist+= temp2*temp2;
		}
		dist = std::sqrt(dist);
	}
}

void interact(simulation* sim){
	(*sim).energy_potential = 0;
	for (int i = 0; i < (*sim).n_types; ++i)
	{
		for (int j = 0; j < (*sim).n_types; ++j)
			{
				interaction[i][j](sim,i,j);
			}		
	}
}

void free_particles(simulation* sim, int type1, int type2)
{
	#pragma omp target teams distribute for
	for (int j = 0; j < (*sim).n_particles[type1]; ++j)
	{
		#pragma omp parallel for
		for (int k = 0; k < n_dimensions; ++k)
		{
			(*sim).acceleration[type1][j][k] == 0;
		}
	}
	#pragma omp target teams distribute for
	for (int j = 0; j < (*sim).n_particles[type2]; ++j)
	{
		#pragma omp parallel for
		for (int k = 0; k < n_dimensions; ++k)
		{
			(*sim).acceleration[type2][j][k] == 0;
		}
	}
}

void lj_periodic(simulation* sim){
	
}