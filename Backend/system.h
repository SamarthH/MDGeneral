#ifndef SYSTEM_H
#define SYSTEM_H 

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

class system_state
{
public:
	double*** position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
	double* temperature; // This defines the temperatures of the n_types particle sets
};

class input_params
{
public:
	int n_types; //This represents the number of types of particles
	int n_dimensions; //This represents the number of dimensions of the simulation (by default must be 3)
	int* n_particles; //This represents the number of particles of each type (to be allocated to an n_types sized array)
	double timestep;
	double runtime;
	int parallelize; // If 0, parallelize. Else, do not parallelize.
	double* mass; //This represents the mass of each type of particle (to be allocated to an n_types sized array)
	int periodic_boundary; //Use periodic boundary conditions if 1. If 0, use rigid walls.
};

class simulation : public system_state, public input_params
{
public:
	void (**thermostat)(input_params*, system_state*); //This stores thermostats for different particle sets
	void (**interaction)(double*, system_state*); //This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.
	double box_size_limits[n_dimensions]; // We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box.
	simulation(int types, int dimensions, int n_par[types], double m[types], int parallel, int periodic, double time, double run){
		n_types = types;
		n_dimensions = dimensions;
		//Allocating and defining n_particles
		n_particles = new int[n_types];
		for (int i = 0; i < n_ty; ++i)
		{
			n_particles[i] = n_par[i];
		}
		//Done
		//Allocating and defining mass
		mass = new int[n_types];
		//Done
		for (int i = 0; i < n_dimensions; ++i)
		{
			box_size_limits[i] = size[i];
		}
	};
	~simulation();
	
};

#endif
