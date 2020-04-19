#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "interaction.h"

double distance_periodic(System::simulation& sim, int type1, int n1, int type2, int n2)
{
	double dist = 0;
	for (int i = 0; i < sim.n_dimensions; ++i)
	{
		double temp = abs(sim.position[type1][n1][i] - sim.position[type2][n2][i]);
		double temp2 = std::min(sim.box_size_limits[i] - temp,temp);
		dist+= temp2*temp2;
	}
	dist = std::sqrt(dist);
	return dist;
}

void interact(System::simulation& sim){
	sim.energy_potential = 0;

	#pragma omp parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp target teams distribute parallel for collapse(2)
		for (int j = 0; j < sim.n_particles[i]; ++j)
		{
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.acceleration[i][j][k] = 0;
			}
		}
	}


	//This implementation is for symmetric interactions only (which makes the most sense)
	for (int i = 0; i < sim.n_types; ++i)
	{
		for (int j = i; j < sim.n_types; ++j)
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

void lj_periodic(System::simulation& sim, int type1, int type2){
	// As of now, assuming that r_c <= min(box_size)/2

	double epot =0; //Temp storage of potential energy
	#pragma omp target teams distribute parallel for collapse(2)
	for (int i = 0; i < sim.n_particles[type1]; ++i)
	{
		for (int j = 0; j < sim.n_particles[type2]; ++j)
		{
			double r = distance_periodic(sim,type1,i,type2,j);
			if(r < sim.rcut_lj)
			{
				double f = 0;
				double r2 = r*r;
				double r6 = r2*r2*r2;

				double b1 = 4*sim.epsilon_lj*sim.sigma_lj_6/r6;
				double b2 = sim.sigma_lj_6/r6;

				epot+= b1*(b2-1);

				f = 6*b1*(2*b2-1)/r2;

				double fx;
				#pragma omp parallel for
				for (int k = 0; k < sim.n_dimensions; ++k)
				{
					double x = sim.position[type1][i][k] - sim.position[type2][j][k];
					if(x > sim.box_size_limits[k])
					{
						x -= -sim.box_size_limits[k];
					}
					else if(x < -sim.box_size_limits[k])
					{
						x += sim.box_size_limits[k];
					}

					fx = f*x;
					#pragma omp atomic
					sim.acceleration[type1][i][k] += fx/sim.mass[type1];
					#pragma omp atomic
					sim.acceleration[type2][j][k] -= fx/sim.mass[type2];
				}
			}
		}
	}


	if(type1 == type2){
		epot/=2;
	}

	sim.energy_potential+=epot;
}