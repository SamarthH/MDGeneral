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

	class molecule_state
	{
	public:

		std::array<double,3> position_com; ///< Position of COM of the molecule
		std::array<double,3> velocity_com; ///< Velocity of COM of the molecule
		std::array<double,3> acceleration_com; ///< Acceleration of COM of the molecule
		
		std::vector<std::array<double,3>> position_par_world; ///< Position of each atom in molecule w.r.t. World frame
		std::vector<std::array<double,3>> position_par_com; ///< Position of each atom in molecule w.r.t. COM
		std::vector<std::array<double,3>> position_par_com_init; ///< Initial position of each atom in molecule w.r.t. World frame
		
		std::vector<std::array<double,3>> force_par; ///< Forces on each atom in the World frame

		quaternion quatrot; ///< Rotation quaternion of the molecule
		std::array<double,3> angmomentum; ///< Angular momentum of the molecule in World frame
		std::array<double,3> torque; ///< Torque on molecule in world frame
		std::array<std::array<double,3>,3> rotation_matrix; ///< Rotation matrix corresponding to quatrot

		molecule_state() {} // Do nothing
		~molecule_state();

		/// Reserves space for the vectors
		reservespace(int n_atoms);
		
	};

	class system_state
	{
	public:

		std::vector<std::vector<molecule_state>> mol_state; ///< n_types X n_molecules[type] vector of molecule states

		std::vector<double> temperature; /**< This defines the temperatures of the n_types particle sets */
		double energy_total; /**< Defines the total energy at this instant */
		double energy_potential; /**< Defines the total potential energy of interaction at this instant */
		std::vector<double> energy_kinetic; /**< Defines the kinetic energy of each particle type */
		double time; /**< This is the amount of time passed since the beginning of the simulation */
		long int state; ///< The timestep number the system is in now
		int numpartot;///< Total number of particles

		system_state(int n_types, std::vector<int>& n_molecules, std::vector<int>& n_atoms);
	};

	class molecule_const
	{
	public:

		int n_atoms;
		double mass;
		std::array<double,3> pos_init_mol; ///< This represents the initial positions of the atoms in the molecule w.r.t the COM
		std::array<std::array<double,3>,3> inertia_tensor; ///< Inertia tensor in initial position
		std::array<double,3> inv_inertia_tensor; ///< Inverse of Inertia tensor in initial position (This is required to be the three eigenvalues in order. inv_inerta_tensor[i] = I_diagonal^-1[i][i]). This decreases the number of multiplications required drastically.

		molecule_const();
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

	class interaction_method
	{
	public:
		/**
		 *  \brief Stores the constants of interaction for all the interactions between particles.
		 *
		 *  Stores the constants of interaction for all the interactions between particles. For interaction of atoms n1 and n2 of molecules of type1 and type2 the constants are stored in vector interaction_const[type1][type2][n1][n2] \n
		   If the interaction is of the Lennard Jones type, \n
		   constant[0] = $\epsilon$ \n
		   constant[1] = $\sigma$ \n
		   constant[2] = Cutoff radius/distance (r_cut) \n
		   constant[3] = Truncated Potential (etrunc) \n
		   constant[4] = $\sigma^6$ \n
		   constant[5] = Tail Energy (assuming constant distribution outside cutoff radius) \n
		*/
		std::vector<double> constant;
		/// The interaction function. We need to pass to f_interaction(simulation, molecule type1, molecule type2, atom n1, atom n2, constant )
		/// This returns the potential energy of interaction excluding the tail energy contribution.
		void (*f_interaction)(System::simulation&, int,int,int,int, std::vector<double>);

		void doInteraction(System::simulation& sim, int type1, int type2, int n1, int n2)
		{
			f_interaction(sim,type1,type2,n1,n2,constant);
		}

		interaction_method();
	};

	class thermostat_method
	{
	public:
		std::vector<double> constant; ///< This vector stores constants of the thermostats

		/// This thermostat function is called as f_thermostat(System::simulation& sim, int type, std::vector<double> constant)
		void (*f_thermostat)(System::simulation&, int , std::vector<double>)

		void doThermostat(System::simulation& sim, int type)
		{
			f_thermostat(sim,type,constant);
		}
		thermostat_method();
	};

	/*
	class correlation
	{
	public:
		std::vector<std::vector<std::array<double,3>>> velocity_initial; ///< Stores the velocity of the particles at t=0
		std::vector<std::vector<double>> correlation_velocity; ///< Stores the velocity correlation for each particle type at each timestep (the n_types+1 th entry is the correlation over all types) \n The format is correlation_velocity[step_number][particletype]

		
		///This constructor reserves space for the correlation arrays and initial conditions
		correlation(int n_types, int n_dimensions, std::vector<int>& n_molecules, double runtime, double timestep);
		~correlation();
		
	};
	*/

	class simulation : public input_params, public system_state
	{
	public:
		
		std::vector<molecule_const> molconst; ///< This stores the molecule constants for each type of molecule

		std::vector<thermostat_method> thermostat; /**< This stores thermostats for different particle sets */
		
		std::vector<std::vector<std::vector<std::vector<interaction_method>>>> interaction; /**< This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.*/
		
		std::array<double,3> box_size_limits; /**< We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box (allocate to n_dimensions size) */
		
		int total_steps; ///< Total number of steps to be taken
		
		std::vector<int> dof; ///< This stores the number of degrees of freedom for each molecule/particle type.

		simulation(std::string input):input_params(input), system_state(n_types,n_molecules,n_atoms)
		{

			total_steps = (int)(runtime/timestep);

			//Allocating functions
			try{
				molconst.reserve(n_atoms);

				dof.reserve(n_types);
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

			//Defining dof
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