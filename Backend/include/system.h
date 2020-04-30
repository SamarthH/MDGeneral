/** @file */ 
#ifndef SYSTEM_H
#define SYSTEM_H 

#include <iostream>
#include <omp.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "algorithm_constants.h"

namespace System{
	class system_state
	{
	public:
		std::vector<std::vector<std::vector<double>>> position; /**< Vector of n_types X n_particles[of each type] X n_dimensions size storing positions of particles */
		std::vector<std::vector<std::vector<double>>> orientation; /**< Vector of n_types X n_particles[of each type] X n_dimensions size storing orientation of particles */
		std::vector<std::vector<std::vector<double>>> velocity; /**< Vector of n_types X n_particles[of each type] X n_dimensions size storing velocity of particles */
		std::vector<std::vector<std::vector<double>>> acceleration; /**< Vector of n_types X n_particles[of each type] X n_dimensions size storing accelerations of particles */
		std::vector<double> temperature; /**< This defines the temperatures of the n_types particle sets */
		double energy_total; /**< Defines the total energy at this instant */
		double energy_potential; /**< Defines the total potential energy of interaction at this instant */
		std::vector<double> energy_kinetic; /**< Defines the kinetic energy of each particle type */
		double time; /**< This is the amount of time passed since the beginning of the simulation */
		int state; ///< The timestep number the system is in now
		int numpartot;///< Total number of particles
		system_state(int n_types, int n_dimensions, std::vector<int>& n_particles){
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

					numpartot += n_particles[i];

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
		}
	};

	class input_params
	{
	public:
		int n_types; ///< This represents the number of types of particles
		int n_dimensions; ///< This represents the number of dimensions of the simulation (by default must be 3)
		std::vector<int> n_particles; ///< This represents the number of particles of each type (n_types sized)
		double timestep; ///< This defines the size of each timestep (dt)
		double runtime; ///< This defines the time for which to run the simulation
		int parallelize; ///< If 1, parallelize. Else, do not parallelize.
		std::vector<double> mass; ///< This represents the mass of each type of particle (n_types sized)
		std::vector<double> temperature_required; ///< This is the vector of the temperatures required to be mainted for each particle type by the thermostat.
		int periodic_boundary; ///< Use periodic boundary conditions if 1. If 0, use rigid walls.
		
		input_params(std::string input){
			std::vector<int> input_vector;
			std::stringstream ss(input);
			
			while(ss.good()){			//This packages the input string (comma separated) into the input_vector object
				std::string temp;
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

	class constants_interaction
	{
	public:
		/**
		 *  \brief Stores the constants of interaction for all the interactions between particles.
		 *
		 *  Stores the constants of interaction for all the interactions between particles. For interaction of particles of type1 and type2 the constants are stored in vector interaction_const[type1][type2] \n
		   If the interaction is of the Lennard Jones type, \n
		   interaction_const[i][j][0] = $\epsilon$ \n
		   interaction_const[i][j][1] = $\sigma$ \n
		   interaction_const[i][j][2] = Cutoff radius/distance (r_cut) \n
		   interaction_const[i][j][3] = Truncated Potential (etrunc) \n
		   interaction_const[i][j][4] = $\sigma^6$ \n
		   interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius) \n
		*/
		std::vector<std::vector<std::vector<double>>> interaction_const;


		// Parametrized Constructor

		constants_interaction(int n_types /** Number of types of particles */)
		{
			try{
				interaction_const.reserve(n_types);
				for (int i = 0; i < n_types; ++i)
				{
					interaction_const[i].reserve(n_types);
					for (int j = 0; j < n_types; ++j)
					{
						interaction_const[i][j].reserve(8);
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
		}
	};

	class constants_thermostat
	{
	public:
		std::vector<std::vector<double>> thermostat_const; ///< This vector stores constants of the thermostats

		// Parametrized Constructor

		constants_thermostat(int n_types)
		{
			try{
				thermostat_const.reserve(n_types);
				for (int i = 0; i < n_types; ++i)
				{
					thermostat_const[i].reserve(4);
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
		}
	};

	class correlation
	{
	public:
		std::vector<std::vector<std::vector<double>>> velocity_initial; ///< Stores the velocity of the particles at t=0
		std::vector<std::vector<double>> correlation_velocity; /**< Stores the velocity correlation for each particle type at each timestep (the n_types+1 th entry is the correlation over all types) \n The format is correlation_velocity[step_number][particletype]*/

		/**********************************************
		 * This constructor reserves space for the correlation arrays and initial conditions
		 */
		correlation(int n_types, int n_dimensions, std::vector<int>& n_particles, double runtime, double timestep)
		{
			// Reserving for initial arrays

			try
			{
				velocity_initial.reserve(n_types);

				for (int i = 0; i < n_types; ++i)
				{
					velocity_initial[i].reserve(n_particles[i]);

					for (int j = 0; j < n_particles[i]; ++j)
					{
						velocity_initial[i][j].reserve(n_dimensions);
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
			// Done


			int n_steps = (int)(runtime/timestep);
			//Reserving for correlations

			try
			{
				correlation_velocity.reserve(n_steps);

				for (int i = 0; i < n_steps; ++i)
				{
					correlation_velocity[i].reserve(n_types+1);
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
		}
		~correlation();
		
	};

	class simulation : public input_params, public system_state, public constants_interaction, public constants_thermostat, public correlation
	{
	public:
		std::vector<void (*)(simulation&, int)> thermostat; /**< This stores thermostats for different particle sets */
		std::vector<std::vector<void (*)(simulation&, int, int)>> interaction; /**< This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.*/
		std::vector<double> box_size_limits; /**< We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box (allocate to n_dimensions size) */
		int total_steps; ///< Total number of steps to be taken

		simulation(std::string input, double size[]):input_params(input), system_state(n_types,n_dimensions,n_particles), constants_interaction(n_types), constants_thermostat(n_types), correlation(n_types,n_dimensions,n_particles,runtime,timestep)
		{

			total_steps = (int)(runtime/timestep);

			//Allocating and defining box_size_limits
			try{
				box_size_limits.reserve(n_dimensions);
			}
			catch(const std::length_error& le){
				std::cerr<<"Error 0001"<<std::endl; 
				exit(0001);
			}
			catch(const std::bad_alloc& ba){
				std::cerr<<"Error 0002"<<std::endl;
				exit(0002);
			}
			for (int i = 0; i < n_dimensions; ++i)
			{
				box_size_limits.push_back(size[i]);
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



		~simulation();
		
	};
}
#endif