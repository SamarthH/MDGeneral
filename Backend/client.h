typedef struct system_state
{
	double*** position; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** orientation; //Needs to be allocated to have n_types X n_particles X n_dimensions size
	double*** velocity; //Needs to be allocated to have n_types X n_particles X n_dimensions size
}system_state;

typedef struct input_params
{
	int n_types; //This represents the number of types of particles
	int n_dimensions;
	int n_paticles[n_types]; //This represents the number of particles of each type
	double timestep;
	double runtime;
	int parallelize; // If 0, parallelize. Else, do not parallelize.
	void (*thermostat[n_types])(struct input_params*, struct system_state*)
}input_params;