class system_state
{
public:
	double*** position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions sizes
};

class input_params
{
public:
	int n_types; //This represents the number of types of particles
	int n_dimensions; //This represents the number of dimensions of the simulation (by default must be 3)
	int n_paticles[n_types]; //This represents the number of particles of each type
	double timestep;
	double runtime;
	int parallelize; // If 0, parallelize. Else, do not parallelize.
	double mass[n_types]; //This represents the mass of each type of particle
	void (*thermostat[n_types])(input_params*, system_state*); //This stores thermostats for different particle sets
	void (*interaction[n_types][n_types])(double*, system_state*); //This defines the set of functions for interaction between different particle types. Also allows for non-symmetric interaction.
	int periodic_boundary; //Use periodic boundary conditions if 1. If 0, use rigid walls.
	double temper_init[n_types]; //This defines the initial temperatures of the n_types particle sets
};

class box : public input_params
{
public:
	double box_size_limits[n_dimensions]; // We assume that the initial limits are all (0,0,0,...,0) to whatever the limits define for a box.
	box(double size[n_dimensions]){
		for (int i = 0; i < n_dimensions; ++i)
		{
			box_size_limits[i] = size[i];
		}
	};
	~box();
	
};

