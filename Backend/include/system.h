#ifndef SYSTEM_H
#define SYSTEM_H 

//#include <stdio.h>
#include <iostream>
//#include <stdlib.h>
#include <omp.h>
#include <string>
#include <sstream>
#include <vector>
#include "thermostat.h"
#include "interaction.h"
namespace system{
	class system_state
	{
	public:
		std::vector<std::vector<std::vector<double>>> position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
		std::vector<std::vector<std::vector<double>>> orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
		std::vector<std::vector<std::vector<double>>> velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
		std::vector<std::vector<std::vector<double>>> acceleration; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
		std::vector<double> temperature; // This defines the temperatures of the n_types particle sets
		double energy_total; //Defines the total energy at this instant
		double energy_potential; //Defines the total potential energy of interaction at this instant
		std::vector<double> energy_kinetic; //Defines the kinetic energy of each particle type
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
		std::vector<int> n_particles; //This represents the number of particles of each type (n_types sized)
		double timestep;
		double runtime;
		extern int parallelize; // If 1, parallelize. Else, do not parallelize.
		std::vector<double> mass; //This represents the mass of each type of particle (n_types sized)
		std::vector<double> temperature_required; //This is the vector of the temperatures required to be mainted for each particle type by the thermostat.
		int periodic_boundary; //Use periodic boundary conditions if 1. If 0, use rigid walls.
		
		input_params(std::string input){
			std::vector<int> input_vector;
			std::stringstream ss(input);
			
			while(ss.good()){			//This packages the input string (comma separated) into the input_vector object
				string temp;
				std::getline(ss, temp, ',');
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
		std::vector<void (*)(simulation&, int)> thermostat; //This stores thermostats for different particle sets
		std::vector<std::vector<void (*)(simulation&, int)>> interaction; //This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.
		std::vector<double> box_size_limits; // We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box (allocate to n_dimensions size)

		simulation(std::string input):input_params(input_params)
		{

			//Please check here @Sweptile

			//Allocating and defining box_size_limits
			try{
				box_size_limits.reserve(n_dimensions);
			}
			catch(const std::length_error& le){
				std::cerr<<"Error 0001"<<std::endl; 
				exit(0001);
			}
			for (int i = 0; i < n_dimensions; ++i)
			{
				box_size_limits.push_back(size[i]);
			}
			//Done

			//Allocating (reserving) system_state variables
			try{
				temperature.reserve(n_types);
				energy_kinetic.reserve(n_types);
				position.reserve(n_types);
				velocity.reserve(n_types);
				orientation.reserve(n_types);
				acceleration.reserve(n_types);

				for (int i = 0; i < n_types; ++i)
				{
					position[i].reserve(n_particles[i]);
					velocity[i].reserve(n_particles[i]);
					orientation[i].reserve(n_particles[i]);
					acceleration[i].reserve(n_particles[i]);

					for (int j = 0; j < n_particles[i]; ++j)
					{
						position[i][j].reserve(n_dimensions);
						velocity[i][j].reserve(n_dimensions);
						orientation[i][j].reserve(n_dimensions);
						acceleration[i][j].reserve(n_dimensions);
					}
				}
			}
			catch(const std::length_error& le){
				std::cerr<<"Error 0001"<<std::endl; 
				exit(0001);
			}
			catch(const std::bad_alloc& ba){
				std::cerr<<"Error 0002"<<std::endl;
				exit(0002);
			}
			//Done

			//Allocating functions
			try{
				thermostat.reserve(n_types);
				interaction.reserve(n_types);
				for (int i = 0; i < n_types; ++i)
				{
					interaction[i].reserve(n_types);
				}
			}
			catch(const std::length_error& le){
				std::cerr<<"Error 0001"<<std::endl; 
				exit(0001);
			}
			catch(const std::bad_alloc& ba){
				std::cerr<<"Error 0002"<<std::endl;
				exit(0002);
			}
			//Done
		};


		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary*/
		/* Temporary functions alloted to thermostat and interaction. Please remove when input is possible*/
		
		for(int i=0; (i< n_types); i++){
			thermostat[i] = no_thermostat;
			for (int j = 0; j < n_types; ++j)
			{
				interaction[i][j] = free_particles;
			}
		}

		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		/* Done*/
		~simulation();
		
	};
}
#endif