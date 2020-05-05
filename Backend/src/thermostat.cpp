#include "thermostat.h"
#include "trng/yarn5.hpp"
#include "trng/normal_dist.hpp"
#include "trng/uniform01_dist.hpp"
#include "trng/gamma_dist.hpp"

void initialize_thermostats(System::simulation& sim){
	for (int i = 0; i < sim.n_types; ++i)
	{
		if(sim.thermostat[i].f_thermostat == no_thermostat){
			// Do nothing
		}
		else if(sim.thermostat[i].f_thermostat == anderson){
			//Precomputing values for anderson thermostat
			sim.thermostat[i].constant[1] = sqrt(BOLTZ_SI*(sim.temperature_required[i])/sim.mass[i]);
		}
		else if(sim.thermostat[i].f_thermostat == bussi){
			//Precomputing values for Bussi thermostat
			double f1 = sim.timestep/sim.thermostat[i].constant[0];
			sim.thermostat[i].constant[1] = 1 + f1 + (f1*f1/2) + (f1*f1*f1/6) + (f1*f1*f1*f1/24);
			sim.thermostat[i].constant[2] = (f1 + (f1*f1/2) + (f1*f1*f1/6) + (f1*f1*f1*f1/24))*std::sqrt(sim.temperature_required[i]/2);
			sim.thermostat[i].constant[3] = std::sqrt(sim.thermostat[i].constant[1]*sim.thermostat[i].constant[2]);
		}
	}
}

void call_thermostat(System::simulation& sim){
	#pragma omp parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.thermostat[i].doThermostat(sim,i);
	}
}

void no_thermostat(System::simulation& sim, int type){
	//Nothing to be done here
}

void anderson(System::simulation& sim, int type, std::vector<double> constant){
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_molecules[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			trng::yarn5 R;

			R.split(sim.total_steps,sim.state);

			R.split(sim.n_types,type);
			
			trng::uniform01_dist<> unif;
			trng::normal_dist<> norm(0,constant[1]);

			int size = omp_get_num_threads();
			int rank = omp_get_thread_num();

			R.split(size,rank);

			if(unif(R) <= sim.thermostat_const[j][0]*(sim.timestep)){
				sim.velocity_com[type][j][k] = norm(R);
			}
		}
	}
}

void bussi(System::simulation& sim,int type, std::vector<double> constant){

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

 	alpha2 = constant[1] + constant[2]*r2sum/sim.energy_kinetic[type] + constant[3]*r[0]/std::sqrt(sim.energy_kinetic[type]);

 	alpha = std::sqrt(alpha2);
	#pragma omp target teams distribute parallel for collapse(2)
	for (int j = 0; j < sim.n_molecules[type]; ++j)
	{
		for (int k = 0; k < sim.n_dimensions; ++k)
		{
			sim.velocity_com[type][j][k] *= alpha;
		}
	}
}