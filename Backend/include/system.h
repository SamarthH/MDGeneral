/** @file */ 
#ifndef SYSTEM_H
#define SYSTEM_H 

#include <iostream>
#include <omp.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <array>
#include "algorithm_constants.h"
#include "quaternion.h"

namespace System{
	class system_state
	{
	public:

		std::vector<std::vector<std::array<double,3>>> position_com; /**< Vector of n_types X n_molecules[of each type] X n_dimensions (<=3) size storing positions of molecule COM */
		std::vector<std::vector<std::array<double,3>>> velocity_com; /**< Vector of n_types X n_molecules[of each type] X n_dimensions (<=3) size storing velocity of molecule COM */
		std::vector<std::vector<std::array<double,3>>> acceleration_com; /**< Vector of n_types X n_molecules[of each type] X n_dimensions (<=3) size storing accelerations of molecule COM */

		std::vector<std::vector<std::vector<std::array<double, 3>>>> position_par_world; ///< Position of each particle w.r.t. the world frame
		std::vector<std::vector<std::vector<std::array<double, 3>>>> position_par_com; ///< Position of each particle w.r.t. the COM frame
		std::vector<std::vector<std::vector<std::array<double, 3>>>> position_par_com_init; ///< Position of each particle w.r.t. the COM frame at t=0. This is constant throughout the simulation.

		std::vector<std::vector<std::vector<std::array<double, 3>>>> force_par; ///< Force on each particle

		std::vector<std::vector<std::array<double,3>>> angvelocity;///< Vector of angular velocities of molecules
		std::vector<std::vector<std::array<double,3>>> angmomentum; ///< Vector of angular momenta of molecules
		std::vector<std::vector<std::array<double,3>>> torque; ///< Vector of torques of molecules
		std::vector<std::vector<std::array<double,4>>> quatrot; ///< Rotation quaternion of molecules
		std::vector<std::vector<std::array<std::array<double,3>,3>>> rotation_matrix; ///< Rotation Matrix generated from the quaternion

		std::vector<double> temperature; /**< This defines the temperatures of the n_types particle sets */
		double energy_total; /**< Defines the total energy at this instant */
		double energy_potential; /**< Defines the total potential energy of interaction at this instant */
		std::vector<double> energy_kinetic; /**< Defines the kinetic energy of each particle type */
		double time; /**< This is the amount of time passed since the beginning of the simulation */
		long int state; ///< The timestep number the system is in now
		int numpartot;///< Total number of particles

		system_state(int n_types, int n_dimensions, std::vector<int>& n_molecules);
	};

	class molecule_const
	{
	public:

		std::vector<int> n_atoms;
		std::vector<std::array<double,3>> pos_init_mol; ///< This represents the initial positions of the atoms in the molecule w.r.t the COM
		std::vector<std::array<std::array<double,3>,3>> inertia_tensor; ///< Inertia tensor in initial position
		std::vector<std::array<double,3>> inv_inertia_tensor; ///< Inverse of Inertia tensor in initial position (This is required to be the three eigenvalues in order. inv_inerta_tensor[i] = I_diagonal^-1[i][i]). This decreases the number of multiplications required drastically.

		molecule_const(int n_types);
		~molecule_const();
	};
	
	class input_params
	{
	public:
		int n_types; ///< This represents the number of types of particles
		int n_dimensions; ///< This represents the number of dimensions of the simulation (by default must be 3)
		std::vector<int> n_molecules; ///< This represents the number of particles of each type (n_types sized)
		double timestep; ///< This defines the size of each timestep (dt)
		double runtime; ///< This defines the time for which to run the simulation
		int parallelize; ///< If 1, parallelize. Else, do not parallelize.
		std::vector<double> mass; ///< This represents the mass of each type of particle (n_types sized)
		std::vector<double> temperature_required; ///< This is the vector of the temperatures required to be mainted for each particle type by the thermostat.
		int periodic_boundary; ///< Use periodic boundary conditions if 1. If 0, use rigid walls.
		
		input_params(std::string input);
	};

	class constants_interaction
	{
	public:
		/**
		 *  \brief Stores the constants of interaction for all the interactions between particles.
		 *
		 *  Stores the constants of interaction for all the interactions between particles. For interaction of atoms n1 and n2 of molecules of type1 and type2 the constants are stored in vector interaction_const[type1][type2][n1][n2] \n
		   If the interaction is of the Lennard Jones type, \n
		   interaction_const[i][j][k][l][0] = $\epsilon$ \n
		   interaction_const[i][j][k][l][1] = $\sigma$ \n
		   interaction_const[i][j][k][l][2] = Cutoff radius/distance (r_cut) \n
		   interaction_const[i][j][k][l][3] = Truncated Potential (etrunc) \n
		   interaction_const[i][j][k][l][4] = $\sigma^6$ \n
		   interaction_const[i][j][k][l][5] = Tail Energy (assuming constant distribution outside cutoff radius) \n
		*/
		std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> interaction_const;


		// Parametrized Constructor

		constants_interaction(int n_types /** Number of types of particles */);
	};

	class constants_thermostat
	{
	public:
		std::vector<std::vector<double>> thermostat_const; ///< This vector stores constants of the thermostats

		// Parametrized Constructor

		constants_thermostat(int n_types);
	};

	class correlation
	{
	public:
		std::vector<std::vector<std::array<double,3>>> velocity_initial; ///< Stores the velocity of the particles at t=0
		std::vector<std::vector<double>> correlation_velocity; /**< Stores the velocity correlation for each particle type at each timestep (the n_types+1 th entry is the correlation over all types) \n The format is correlation_velocity[step_number][particletype]*/

		/**********************************************
		 * This constructor reserves space for the correlation arrays and initial conditions
		 */
		correlation(int n_types, int n_dimensions, std::vector<int>& n_molecules, double runtime, double timestep);
		~correlation();
		
	};

	class mixingclass1 : public constants_interaction, public constants_thermostat, public molecule_const
	{
	public:
		mixingclass1(int n_types): constants_interaction(n_types), constants_thermostat(n_types), molecule_const(n_types) {}
		~mixingclass1();
		
	};

	class mixingclass2 : public mixingclass1, public system_state, public correlation
	{
	public:
		mixingclass2(int n_types,int n_dimensions,std::vector<int>& n_molecules, double runtime, double timestep) : mixingclass1(n_types), system_state(n_types,n_dimensions,n_molecules), correlation(n_types,n_dimensions,n_molecules,runtime,timestep) {}
		~mixingclass2();
		
	};



	class simulation : public input_params, public mixingclass2
	{
	public:
		std::vector<void (*)(simulation&, int)> thermostat; /**< This stores thermostats for different particle sets */
		std::vector<std::vector<std::vector<std::vector<void (*)(simulation&, int, int, int, int)>>>> interaction; /**< This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.*/
		std::vector<double> box_size_limits; /**< We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box (allocate to n_dimensions size) */
		int total_steps; ///< Total number of steps to be taken
		std::vector<int> dof; ///< This stores the number of degrees of freedom for each molecule/particle type.

		simulation(std::string input, double size[]):input_params(input), mixingclass2(n_types,n_dimensions,n_molecules,runtime,timestep)
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

			//Allocating and defining dof
			try{
				dof.reserve(n_types);
			}
			catch(const std::length_error& le){
				std::cerr<<"Error 0001"<<std::endl; 
				exit(0001);
			}
			catch(const std::bad_alloc& ba){
				std::cerr<<"Error 0002"<<std::endl;
				exit(0002);
			}
			for (int i = 0; i < n_types; ++i)
			{
				dof[i] = n_dimensions;
			}
			//Done
		};



		~simulation();
		
	};
}
#endif