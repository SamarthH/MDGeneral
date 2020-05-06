#include "system.h"
#include "interaction.h"
#include "thermostat.h"
#include "universal_functions.h"

System::system_state::system_state(int n_types, std::vector<int>& n_molecules, std::vector<int>& n_atoms)
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


/*
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
*/
/*
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
*/

void System::molecule_state::reservespace(int n_atoms)
{
	try
	{
		position_par_world.reserve(n_atoms);
		position_par_com.reserve(n_atoms);
		position_par_com_init.reserve(n_atoms);
		force_par.reserve(n_atoms);
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

System::interaction_list::interaction_list()
{
	try
	{
		interactions.reserve(4);
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

void System::interaction_list::doInteractions(System::simulation& sim, int type1, int type2, int n1, int n2)
{
	for (int i = 0; i < n_interactions; ++i)
	{
		interactions[i].doInteraction(sim,type1,type2,n1,n2);
	}
}

void System::interaction_list::getNumberOfInteractions()
{
	n_interactions = interactions.size();
}

void System::interaction_list::initializeConstantArrays(System::simulation& sim, int i, int j)
{
	for (int k = 0; k < n_interactions; ++k)
	{
		interactions[k].initialize_constant_array(sim,i,j);
	}
}

System::interaction_method::interaction_method()
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

void System::interaction_method::doInteraction(System::simulation& sim, int type1, int type2, int n1, int n2)
{
	f_interaction(sim,type1,type2,n1,n2,constant);
}

void System::interaction_method::initialize_constant_array(System::simulation& sim, int i, int j)
{
	int dim = sim.n_dimensions;

	if(f_interaction == lj_periodic || f_interaction == lj_box)
	{
		/*
		constant[0] = \epsilon
		constant[1] = \sigma
		constant[2] = Cutoff radius/distance (r_cut)
		constant[3] = Truncated Potential (etrunc)
		constant[4] = \sigma^6
		constant[5] = Tail Energy (assuming constant distribution outside cutoff radius)
		*/
		constant[4] = std::pow(constant[1],6);
		double temp = constant[4]/std::pow(constant[2],6); //temp = (sigma/r_cut)^6
		constant[3] = 4*constant[0]*temp*(temp-1);
		constant[5] = 2*constant[0]*(sim.n_molecules[i]*sim.n_molecules[j]/sim.vol)*surface_unit_sphere(dim);
		constant[5] *= (constant[4]*constant[4]*std::pow(constant[2],dim-12)/(dim-12) - constant[4]*std::pow(constant[2],dim-6)/(dim-6));

	}
	else
	{
		std::cerr<<"Error 0003"<<std::endl;
		exit(0003);
	}
}