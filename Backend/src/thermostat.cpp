#include "thermostat.h"

void call_thermostat(System::simulation& sim){
	#pragma omp parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.thermostat[i](sim,i);
	}
}

void no_thermostat(System::simulation& sim, int type){
	//Nothing to be done here
}

void anderson(System::simulation& sim, int type){
	std:: random_device rd;
	std::mt19937 gen(rd());

	double stdev_boltzman = sqrt(BOLTZ_SI*(sim.temperature_required[type])/sim.mass[type]);

	std::uniform_real_distribution<> unif(0, 1);
	std::normal_distribution<> norm(0, stdev_boltzman);

	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			/*
			if(unif(gen) <= ANDERSON_NU*(sim.timestep)){
				sim.velocity[type][j][k] = norm(gen);
			}
			*/
			//This has been removed until thread safe alternatives are found and implemented.
		}
	}
}