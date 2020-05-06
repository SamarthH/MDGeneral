/** @file */ 
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "universal_functions.h"
#include "interaction.h"
#include "quaternion.h"

//
//
//
// MAINTAIN DISTANCE
//
//
//
/* Internal Functions start here*/

//Updates the acceleration array
void _get_acceleration(System::simulation& sim)
{
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{

			double f[3] = {0,0,0}; // Force array

			#pragma omp parallel for
			for (int k = 0; k < sim.n_atoms[i]; ++k)
			{
				#pragma omp parallel for
				for (int l = 0; l < sim.n_dimensions; ++l)
				{
					#pragma omp atomic
					f[l] += sim.mol_state[i][j].force_par[k][l];
				}
			}

			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.mol_state[i][j].acceleration_com[k] = f[k]/sim.molconst[i].mass;
			}

		}
	}
}

//Updates the torque array
void _get_torque(System::simulation& sim)
{
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_atoms[i]; ++k)
			{
				std::array<double,3> t;
				qua_vcross(sim.mol_state[i][j].position_par_com[k],sim.mol_state[i][j].force_par[k],t);
				#pragma omp parallel for
				for (int l = 0; l < 3; ++l)
				{
					#pragma omp atomic
					sim.mol_state[i][j].torque[l] += t[l];
				}				
			}
		}
	}
}
/* Internal Functions end here*/
//
//
//
//
// MAINTAIN DISTANCE
//
//
//
//
//


void initialize_interactions(System::simulation& sim)
{
	int dim = sim.n_dimensions;
	double vol = 1;
	for (int i = 0; i < sim.n_dimensions; ++i)
	{
		vol*= sim.box_size_limits[i];
	}
	for (int i = 0; i < sim.n_types; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			for (int k = 0; k < sim.n_atoms[i]; ++k)
			{
				for (int l = 0; l < sim.n_atoms[j]; ++l)
				{
					sim.interaction[i][j][k][l].initializeConstantArrays(sim,i,j);
				}
			}
		}
	}
}

double distance_periodic(System::simulation& sim, int type1, int n1, int m1, int type2, int n2, int m2)
{
	double dist = 0;
	for (int i = 0; i < sim.n_dimensions; ++i)
	{
		double temp = abs(sim.mol_state[type1][n1].position_par_world[m1][i] - sim.mol_state[type2][n2].position_par_world[m2][i]);
		double temp2 = std::min(sim.box_size_limits[i] - temp,temp);
		dist+= temp2*temp2;
	}
	dist = std::sqrt(dist);
	return dist;
}

void interact(System::simulation& sim){
	sim.energy_potential = 0;

	//This implementation is for symmetric interactions only (which makes the most sense)
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j <= i; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_atoms[i]; ++k)
			{
				#pragma omp parallel for
				for (int l = 0; l < sim.n_atoms[j]; ++l)
				{
					sim.interaction[i][j][k][l].doInteractions(sim,i,j,k,l);
				}
			}
		}		
	}

	_get_acceleration(sim);
	_get_torque(sim);
}

void free_particles(System::simulation& sim, int type1, int type2, std::vector<double> constant)
{
	//This does nothing. Don't worry
}

void lj_periodic(System::simulation& sim, int type1, int type2, int k, int l, std::vector<double> constant){
	
	double epot =0; //Temp storage of potential energy
	#pragma omp target teams distribute parallel for collapse(2)
	for (int i = 0; i < sim.n_molecules[type1]; ++i)
	{
		for (int j = 0; j < sim.n_molecules[type2]; ++j)
		{
			double r = distance_periodic(sim,type1,i,k,type2,j,l);
			if(r < constant[2])
			{
				double f = 0;
				double r2 = r*r;
				double r6 = r2*r2*r2;

				double b1 = 4*constant[0]*constant[4]/r6;
				double b2 = constant[4]/r6;

				f = 6*b1*(2*b2-1)/r2;

				#pragma omp parallel for
				for (int m = 0; m < sim.n_dimensions; ++m)
				{
					double fx;
					double x = sim.mol_state[type1][i].position_par_world[k][m] - sim.mol_state[type2][j].position_par_world[l][m];
					if(x > sim.box_size_limits[m])
					{
						x -= -sim.box_size_limits[m];
					}
					else if(x < -sim.box_size_limits[m])
					{
						x += sim.box_size_limits[m];
					}

					fx = f*x;
					#pragma omp atomic
					sim.mol_state[type1][i].force_par[k][m] += fx;
					#pragma omp atomic
					sim.mol_state[type2][j].force_par[l][m] -= fx;
				}
				#pragma omp atomic
				epot+= (b1*(b2-1)-constant[3]);
			}
		}
	}


	if(type1 == type2){
		epot/=2;
	}

	epot+=constant[5];

	#pragma omp atomic
	sim.energy_potential+=epot;
}

void lj_box(System::simulation& sim, int type1, int type2, int k, int l, std::vector<double> constant){

	/*
	interaction_const[i][j][0] = \epsilon
	interaction_const[i][j][1] = \sigma
	interaction_const[i][j][2] = Cutoff radius/distance (r_cut)
	interaction_const[i][j][3] = Truncated Potential (etrunc)
	interaction_const[i][j][4] = \sigma^6
	interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius)
	*/
	double epot =0; //Temp storage of potential energy
	#pragma omp target teams distribute parallel for collapse(2)
	for (int i = 0; i < sim.n_molecules[type1]; ++i)
	{
		for (int j = 0; j < sim.n_molecules[type2]; ++j)
		{
			double r2 = 0;
			for (int m = 0; m < sim.n_dimensions; ++m)
			{
				r2 += std::pow(sim.mol_state[type1][i].position_par_world[k][m] - sim.mol_state[type2][j].position_par_world[l][m],2);
			}
			if(std::sqrt(r2) < constant[2])
			{
				double f = 0;
				double r6 = r2*r2*r2;

				double b1 = 4*constant[0]*constant[4]/r6;
				double b2 = constant[4]/r6;

				f = 6*b1*(2*b2-1)/r2;

				#pragma omp parallel for
				for (int m = 0; m < sim.n_dimensions; ++m)
				{
					double fx;
					double x = sim.mol_state[type1][i].position_par_world[k][m] - sim.mol_state[type2][j].position_par_world[l][m];
					fx = f*x;
					#pragma omp atomic
					sim.mol_state[type1][i].force_par[k][m] += fx;
					#pragma omp atomic
					sim.mol_state[type2][j].force_par[l][m] -= fx;
				}
				#pragma omp atomic
				epot+= (b1*(b2-1)-constant[3]);
			}
		}
	}


	if(type1 == type2){
		epot/=2;
	}

	epot+=constant[5];

	#pragma omp atomic
	sim.energy_potential+=epot;
}