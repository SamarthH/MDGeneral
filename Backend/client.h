class system_state
{
	double*** position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions size
};

class input_params
{
	int n_types; //This represents the number of types of particles
	int n_dimensions;
	int n_paticles[n_types]; //This represents the number of particles of each type
	double timestep;
	double runtime;
	int parallelize; // If 0, parallelize. Else, do not parallelize.
	double mass[n_types]; //This represents the mass of each type of particle
	void (*thermostat[n_types])(input_params*, system_state*); //This stores thermostats for different particle sets
	void (*interaction[n_types][n_types])(double*, system_state*);
};

