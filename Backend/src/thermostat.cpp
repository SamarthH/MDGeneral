#include "thermostat.h"
#include "trng/yarn2.hpp"
#include "trng/normal_dist.hpp"
#include "trng/uniform01_dist.hpp"

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

	double stdev_boltzman = sqrt(BOLTZ_SI*(sim.temperature_required[type])/sim.mass[type]);

	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			trng::yarn2 R1,R2;

			R1.split(sim.n_types,type);
			R2.split(sim.n_types,type);

			R1.split(2,0);
			R2.split(2,1);
			
			trng::uniform01_dist<> unif;
			trng::normal_dist<> norm(0,stdev_boltzman);

			int size = omp_get_num_threads();
			int rank = omp_get_thread_num();

			R1.split(size,rank);
			R2.split(size,rank);

			if(unif(R1) <= sim.anderson_nu*(sim.timestep)){
				sim.velocity[type][j][k] = norm(R2);
			}
		}
	}
}