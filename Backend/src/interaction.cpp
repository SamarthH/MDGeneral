/** @file */ 
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "universal_functions.h"
#include "interaction.h"

/*******************************************************************************
 * \brief Initializes the constant arrays for interactions for speed
 * 
 * The function initializes the arrays for more efficient computation by precomputing the required factors 
 *
 * @param sim Simulation being initialized
 ******************************************************************************/
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
		for (int j = 0; j < sim.n_types; ++j)
		{
			if(sim.interaction[i][j] == lj_periodic || sim.interaction[i][j] == lj_box)
			{
				/*
				interaction_const[i][j][0] = \epsilon
				interaction_const[i][j][1] = \sigma
				interaction_const[i][j][2] = Cutoff radius/distance (r_cut)
				interaction_const[i][j][3] = Truncated Potential (etrunc)
				interaction_const[i][j][4] = \sigma^6
				interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius)
				*/
				sim.interaction_const[i][j][4] = std::pow(sim.interaction_const[i][j][1],6);
				double temp = sim.interaction_const[i][j][4]/std::pow(sim.interaction_const[i][j][2],6); //temp = (sigma/r_cut)^6
				sim.interaction_const[i][j][3] = 4*sim.interaction_const[i][j][0]*temp*(temp-1);
				sim.interaction_const[i][j][5] = 2*sim.interaction_const[i][j][0]*(sim.n_particles[i]*sim.n_particles[j]/vol)*surface_unit_sphere(dim);
				sim.interaction_const[i][j][5] *= (sim.interaction_const[i][j][4]*sim.interaction_const[i][j][4]*std::pow(sim.interaction_const[i][j][2],dim-12)/(dim-12) - sim.interaction_const[i][j][4]*std::pow(sim.interaction_const[i][j][2],dim-6)/(dim-6));

			}
			else
			{
				std::cerr<<"Error 0003"<<std::endl;
				exit(0003);
			}
		}
	}
}

/*******************************************************************************
 * \brief Returns the distance between two particles for periodic boundary conditions
 * 
 * Returns the distance between two particles labelled number n1,n2 of type1,type2 
 * respectively in the position array assuming periodic boundary conditions
 *
 * @param sim Simulation being used
 * @param type1 Particle type of particle 1
 * @param type2 Particle type of particle 2
 * @param n1 Particle index of particle 1 in position[type1] array
 * @param n2 Particle index of particle 2 in position[type2] array
 ******************************************************************************/
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

/*******************************************************************************
 * \brief This function calls all the required interaction functions between the particles
 * 
 * Calls all the interactiosn and updates acceleration and energy arrays
 *
 * @param sim Simulation being used
 ******************************************************************************/
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

/*******************************************************************************
 * \brief Setup the free particle interaction between two particle types
 * 
 * Setup of free particle interaction between particle types type1 and type2.
 * This does nothing. The function is empty.
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void free_particles(System::simulation& sim, int type1, int type2)
{
	//This does nothing. Don't worry
}

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for periodic boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for periodic boundary conditions between particle types type1 and type2.
 * This requires the cutoff <= half the box size because it only checks the nearest images.
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void lj_periodic(System::simulation& sim, int type1, int type2){
	// As of now, assuming that r_c <= min(box_size)/2

	/*
	interaction_const[i][j][0] = \epsilon
	interaction_const[i][j][1] = \sigma
	interaction_const[i][j][2] = Cutoff radius/distance (r_cut)
	interaction_const[i][j][3] = Truncated Potential (etrunc)
	interaction_const[i][j][4] = \sigma^6
	interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius)
	*/


	double epot =0; //Temp storage of potential energy
	#pragma omp target teams distribute parallel for collapse(2) reduction(+ : epot)
	for (int i = 0; i < sim.n_particles[type1]; ++i)
	{
		for (int j = 0; j < sim.n_particles[type2]; ++j)
		{
			double r = distance_periodic(sim,type1,i,type2,j);
			if(r < sim.interaction_const[type1][type2][2])
			{
				double f = 0;
				double r2 = r*r;
				double r6 = r2*r2*r2;

				double b1 = 4*sim.interaction_const[type1][type2][0]*sim.interaction_const[type1][type2][4]/r6;
				double b2 = sim.interaction_const[type1][type2][4]/r6;

				epot+= (b1*(b2-1)-sim.interaction_const[type1][type2][3]);

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

	epot+=sim.interaction_const[type1][type2][5];

	#pragma omp atomic
	sim.energy_potential+=epot;
}

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for rigid box boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for rigid box boundary conditions between particle types type1 and type2.
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void lj_box(System::simulation& sim, int type1, int type2){

	/*
	interaction_const[i][j][0] = \epsilon
	interaction_const[i][j][1] = \sigma
	interaction_const[i][j][2] = Cutoff radius/distance (r_cut)
	interaction_const[i][j][3] = Truncated Potential (etrunc)
	interaction_const[i][j][4] = \sigma^6
	interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius)
	*/

	double epot =0; //Temp storage of potential energy
	#pragma omp target teams distribute parallel for collapse(2) reduction(+ : epot)
	for (int i = 0; i < sim.n_particles[type1]; ++i)
	{
		for (int j = 0; j < sim.n_particles[type2]; ++j)
		{
			double r2 = 0;

			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				r2 += std::pow(sim.position[type1][i][k] - sim.position[type2][j][k],2);
			}

			if(std::sqrt(r2) < sim.interaction_const[type1][type2][2])
			{
				double f = 0;
				double r6 = r2*r2*r2;

				double b1 = 4*sim.interaction_const[type1][type2][0]*sim.interaction_const[type1][type2][4]/r6;
				double b2 = sim.interaction_const[type1][type2][4]/r6;

				epot+= (b1*(b2-1)-sim.interaction_const[type1][type2][3]);

				f = 6*b1*(2*b2-1)/r2;

				double fx;
				#pragma omp parallel for
				for (int k = 0; k < sim.n_dimensions; ++k)
				{
					double x = sim.position[type1][i][k] - sim.position[type2][j][k];
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

	epot+=sim.interaction_const[type1][type2][5];

	#pragma omp atomic
	sim.energy_potential+=epot;
}