#include "thermostat.h"
#include "trng/yarn5.hpp"
#include "trng/normal_dist.hpp"
#include "trng/uniform01_dist.hpp"
#include "trng/gamma_dist.hpp"

void initialize_thermostats(System::simulation& sim){
	for (int i = 0; i < n_types; ++i)
	{
		if(sim.thermostat[i] == no_thermostat){
			// Do nothing
		}
		else if(sim.thermostat[i] == anderson){
			//Precomputing values for anderson thermostat
			sim.thermostat_const[i][1] = sqrt(BOLTZ_SI*(sim.temperature_required[i])/sim.mass[i]);
		}
		else if(sim.thermostat[i] == bussi){
			//Precomputing values for Bussi thermostat
			double f1 = sim.timestep/sim.thermostat_const[i][0];
			sim.thermostat_const[i][1] = 1 + f1 + (f1*f1/2) + (f1*f1*f1/6) + (f1*f1*f1*f1/24);
			sim.thermostat_const[i][2] = (f1 + (f1*f1/2) + (f1*f1*f1/6) + (f1*f1*f1*f1/24))*std::sqrt(sim.temperature_required[i]/2);
			sim.thermostat_const[i][3] = std::sqrt(a*b);
		}
	}
}

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
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			trng::yarn5 R;

			R.split(sim.total_steps,sim.state);

			R.split(sim.n_types,type);
			
			trng::uniform01_dist<> unif;
			trng::normal_dist<> norm(0,sim.thermostat_const[i][1]);

			int size = omp_get_num_threads();
			int rank = omp_get_thread_num();

			R.split(size,rank);

			if(unif(R) <= sim.thermostat_const[j][0]*(sim.timestep)){
				sim.velocity[type][j][k] = norm(R);
			}
		}
	}
}

void bussi(System::simulation& sim,int type){

	trng::yarn5 rang;

	rang.split(sim.total_steps,sim.state);

	rang.split(sim.n_types,type);
	
	trng::normal_dist<> norm(0,1);

	double alpha2,alpha;

	double r[sim.dof[type]];

	for (int i = 0; i < sim.dof[type]; ++i)
	{
		r[i] = norm(rang);
 	}

 	double r2sum = 0;
 	for (int i = 0; i < sim.dof[type]; ++i)
 	{
 		r2sum += r[i]*r[i];
 	}

 	alpha2 = thermostat_const[type][1] + thermostat_const[type][2]*r2sum/sim.energy_kinetic[type] + thermostat_const[type][3]*r[0]/std::sqrt(sim.energy_kinetic[type]);

 	alpha = std::sqrt(alpha2);
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_particles[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			sim.velocity[type][k] *= alpha;
		}
	}
}