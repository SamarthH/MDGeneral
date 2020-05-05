#include "system.h"

System::system_state::system_state(int n_types, int n_dimensions, std::vector<int>& n_molecules)
{
	//Allocating (reserving) system_state variables
	try{
		temperature.reserve(n_types);
		energy_kinetic.reserve(n_types);
		mol_state.reserve(n_types);

		for (int i = 0; i < n_types; ++i)
		{
			mol_state[i].reserve(n_molecules[i]);

			for (int j = 0; j < n_molecules[i]; ++j)
			{
				mol_state[i][j].reservespace(n_atoms[i]);
			}

			numpartot += n_molecules[i];
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

System::input_params::input_params(std::string input)
{
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
	n_molecules.resize(n_types);
	for(int i=2; i<n_types+2; i++){
		n_molecules.push_back((int)input_vector[i]);
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

System::constants_interaction::constants_interaction()
{
	try{
		constant.reserve(8);
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

System::constants_thermostat::constants_thermostat(int n_types)
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

System::correlation::correlation(int n_types, int n_dimensions, std::vector<int>& n_molecules, double runtime, double timestep)
{
	// Reserving for initial arrays

	try
	{
		velocity_initial.reserve(n_types);

		for (int i = 0; i < n_types; ++i)
		{
			velocity_initial[i].reserve(n_molecules[i]);
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

System::molecule_const::molecule_const(int n_types)
{
	try{
		n_atoms.reserve(n_types);
		pos_init_mol.reserve(n_types);
		inertia_tensor.reserve(n_types);
		inv_inertia_tensor.reserve(n_types);
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