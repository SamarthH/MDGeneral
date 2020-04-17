#ifndef SYSTEM_H
#define SYSTEM_H 

//#include <stdio.h>
#include <iostream>
//#include <stdlib.h>
#include <omp.h>
#include <string>
#include <sstream>
#include <vector>

class system_state
{
public:
	double*** position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
	double*** acceleration; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
	double* temperature; // This defines the temperatures of the n_types particle sets
	double energy_total; //Defines the total energy at this instant
	double energy_potential; //Defines the total potential energy of interaction at this instant
	double* energy_kinetic; //Defines the kinetic energy of each particle type
	double time;
	int state;
	system_state(){
	}
};

class input_params
{
public:
	int n_types; //This represents the number of types of particles
	int n_dimensions; //This represents the number of dimensions of the simulation (by default must be 3)
	//int* n_particles; //This represents the number of particles of each type (to be allocated to an n_types sized array)
	std::vector<int> n_particles; //This represents the number of particles of each type (n_types sized)
	double timestep;
	double runtime;
	extern int parallelize; // If 1, parallelize. Else, do not parallelize.
	//double* mass; //This represents the mass of each type of particle (to be allocated to an n_types sized array)
	std::vector<double> mass; //This represents the mass of each type of particle (n_types sized)
	//double* temperature_required; //This is the array of temperature required to be mainted for each particle type by the thermostat.
	std::vector<double> temperature_required; //This is the vector of the temperatures required to be mainted for each particle type by the thermostat.
	int periodic_boundary; //Use periodic boundary conditions if 1. If 0, use rigid walls.
	
	input_params(std::string input){
		std:vector<int> input_vector;
		std::stringstream ss(input);
		
		while(ss.good()){			//This packages the input string (comma separated) into the input_vector object
			string temp;
			std::getline(ss, temp, ",");
			input_vector.push_back(std::stod(temp));
		}
		
		//input_vector to the variables
		n_types = (int)input_vector[0];
		n_dimensions = (int)input_vector[1];
		n_particles.resize(n_types);
		for(int i=2; i<n_types+2; i++){
			n_particles.push_back((int)input_vector[i]);
		}
		timestep = input_vector[n_types+2];
		runtime = input_vector[n_types+3];
		parallelize = (int)input_vector[n_types+4]; 
		mass.resize(n_types);
		for(int i=n_types+5; i<2*n_types+5; i++){
			mass.push_back(input_vector[i]);
		}
		temperature_required.resize(n_types);
		for(int i=2*n_types+5; i<3*n_types+5; i++){
			temperature_required.push_back(input_vector[i]);
		}
		periodic_boundary = (int)input_vector[3*n_types+5];
		
	}
};

class simulation : public system_state, public input_params
{
public:
	void (**thermostat)(simulation*,int); //This stores thermostats for different particle sets
	void (***interaction)(simulation*,int,int); //This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.
	double* box_size_limits; // We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box (allocate to n_dimensions size)

	simulation(int types, int dimensions, int n_par[types], double m[types], int parallel, int periodic, double time, double run){
		state=0;
		n_types = types;
		n_dimensions = dimensions;
		//Allocating and defining n_particles
		n_particles = new int[n_types];
		if(!n_particles){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
		for (int i = 0; i < n_types; ++i)
		{
			n_particles[i] = n_par[i];
		}
		//Done
		//Allocating and defining mass
		mass = new int[n_types];
		if(!mass){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
		for (int i = 0; i < n_types; ++i)
		{
			mass[i] = m[i];
		}
		//Done
		//Allocating and defining box_size_limits
		box_size_limits = new int[n_dimensions];
		if(!box_size_limits){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
		for (int i = 0; i < n_dimensions; ++i)
		{
			box_size_limits[i] = size[i];
		}
		//Done
		//Allocating space for position, orientation and velocity array (still would need to be randomly innitialized)
		position = new double**[n_types];
		velocity = new double**[n_types];
		orientation = new double**[n_types];
		acceleration = new double**[n_types];
		if(!position || !velocity || !orientation){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
		for (int i = 0; i < n_types; ++i)
		{
			position[i] = new double*[n_particles[i]];
			velocity[i] = new double*[n_particles[i]];
			orientation[i] = new double*[n_particles[i]];
			acceleration[i] = new double*[n_particles[i]];
			if(!position[i] || !velocity[i] || !orientation[i] || !acceleration[i]){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
			for (int j = 0; j < n_particles[i]; ++j)
			{
				position[i][j] = new double*[n_dimensions];
				velocity[i][j] = new double*[n_dimensions];
				orientation[i][j] = new double*[n_dimensions];
				acceleration[i][j] = new double*[n_dimensions];
				if(!position[i][j] || !velocity[i][j] || !orientation[i][j] || !acceleration[i][j]){std::cerr<<"Error 0001"<<std::endl; exit(0001);}
			}
		}

		//Done

		//Allocating and initializing temperature, energy_kinetic

		temperature = new double[n_types];
		energy_kinetic = new double[n_types];

		if(!temperature || !energy_kinetic){std::cerr<<"Error 0001"<<std::endl; exit(0001);}

		for (int i = 0; i < n_types; ++i)
		{
			temperature[i] = energy_kinetic[i] = 0;
		}
		//Done
	};
	~simulation();
	
};

class constants_lj //Only need to call such a class if the interaction to be used is Lennard Jones (LJ) potential. Also, remember that 2.5*\sigma < min(box_size_limits)
{

public:
	constants_lj();
	~constants_lj();
	
};

#endif
