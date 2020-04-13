#include "thermostat.h"

void call_thermostat(simulation* sim){
	#pragma omp parallel for
	for (int i = 0; i < n_type; ++i)
	{
		anderson(sim,i);
	}
}

void anderson(simulation* sim, int type){
	std:: random_device rd;
	std::mt19937 gen(rd());

	double stdev_boltzman = sqrt(BOLTZMAN_SI*((*sim).temperature_required[type])/mass[type]);

	std::uniform_real_distribution<> unif(0, 1);
	std::normal_distribution<> norm(0, stdev_bolzman);

	#pragma omp target teams distribute for
	for (int j = 0; j < (*sim).n_particles[type]; ++j)
	{
		#pragma omp parallel for
		for (int k = 0; k < n_dimensions; ++k)
		{
			if(unif(gen) <= ANDERSON_NU*((*sim).timestep)){
				(*sim).velocity[type][j][k] = norm(gen);
			}
		}
	}
}