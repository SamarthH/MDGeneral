/** @file */ 
#include "integrate.h"
#include "quaternion.h"

//
//
//
// MAINTAIN DISTANCE
//
//
//
/* Internal Functions start here*/

void _get_rotation_matrix(System::simulation& sim)
{
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			generate_rotmat(sim.quatrot[i][j],sim.rotation_matrix[i][j]);
		}
	}
}

void _verlet_trans1(System::simulation& sim)
{
	double dt = sim.timestep;
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.position_com[i][j][k] += dt*sim.velocity_com[i][j][k] + 0.5*dt*dt*sim.acceleration_com[i][j][k];
				sim.position_com[i][j][k] -= sim.box_size_limits[k]*std::floor(sim.position_com[i][j][k]/sim.box_size_limits[k]);
				sim.velocity_com[i][j][k] += dt*sim.acceleration_com[i][j][k]/2;
			}
		}
	}
}

void _verlet_trans2(System::simulation& sim)
{
	double dt = sim.timestep;
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.energy_kinetic[i] = 0;
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.velocity_com[i][j][k] += dt*sim.acceleration_com[i][j][k]/2;
				#pragma omp atomic
				sim.energy_kinetic[i] += sim.velocity_com[i][j][k]*sim.velocity_com[i][j][k];
			}
		}
		sim.energy_kinetic[i] *= sim.mass[i];
		sim.temperature[i] = sim.energy_kinetic[i]/(sim.n_molecules[i]*BOLTZ_SI*sim.n_dimensions);
		sim.energy_total += sim.energy_kinetic[i];
	}
}

/*
 Rotational verlet has been implemented according to section IV A of Robust rotational-velocity-Verlet integration methods, Dmitri Rozmanov and Peter G. Kusalik
 http://dx.doi.org/10.1103/PhysRevE.81.056706
*/
void _verlet_rot1(System::simulation& sim)
{
	double dt = sim.timestep;
	double dt_half = sim.timestep/2;
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			std::array<double,3> L_m, T_m; // These are L and \tau in molecule frame of reference
			std::array<std::array<double,3>,3> Rt;
			transpose_mat(sim.rotation_matrix[i][j],Rt);
			matrix_mult_vec(Rt,sim.angmomentum[i][j],L_m);
			matrix_mult_vec(Rt,sim.torque[i][j],T_m);

			std::array<double,3> omega;
			for (int k = 0; k < 3; ++k)
			{
				omega[k] = sim.inv_inertia_tensor[i][k]*L_m[k];
			}

			std::array<double,3> L_m_dot;
			qua_vcross(omega,L_m,L_m_dot);
			for (int k = 0; k < 3; ++k)
			{
				L_m_dot[k] = T_m[k] - L_m_dot[k];
			}

			for (int k = 0; k < 3; ++k)
			{
				L_m[k] = L_m[k] + (dt_half* L_m_dot[k]); //Updating to L_m(0.5dt)
			}

			quaternion q0, q_half_old, q_half_new;
			for (int k = 0; k < 4; ++k)
			{
				q_half_old[k] = 0;
			}

			//Updating omega
			for (int k = 0; k < 3; ++k)
			{
				omega[k] = sim.inv_inertia_tensor[i][k]*L_m[k];
			}

			quaternion qdot;
			qua_vec_mult(q0,omega,qdot);
			for (int k = 0; k < 4; ++k)
			{
				qdot[k] /= 2;
			}

			for (int k = 0; k < 4; ++k)
			{
				q_half_new[k] = q0[k] + dt_half*qdot[k];
			}

			qua_normalize(q_half_new);

			//Updating the world angular momenta

			for (int k = 0; k < 3; ++k)
			{
				sim.angmomentum[i][j][k] += dt_half*sim.torque[i][j][k];
			}

			quaternion diffq;
			double diff;
			qua_sub(q_half_new,q_half_old,diffq);
			diff = qua_mod(diffq);

			//Starting iterative scheme for finding qdot

			while(diff > QUATERNION_ERROR_TOLERANCE)
			{
				q_half_old = q_half_new; //Copying into old

				qua_rot_vec_opp(q_half_new,sim.angmomentum[i][j],L_m);

				for (int k = 0; k < 3; ++k)
				{
					omega[k] = sim.inv_inertia_tensor[i][k]*L_m[k];
				}

				qua_vec_mult(q_half_new,omega,qdot);
				for (int k = 0; k < 4; ++k)
				{
					qdot[k] /= 2;
				}

				for (int k = 0; k < 4; ++k)
				{
					q_half_new[k] = q0[k] + dt_half*qdot[k];
				}
				qua_normalize(q_half_new);

				qua_sub(q_half_new,q_half_old,diffq);
				diff = qua_mod(diffq);
			}

			for (int k = 0; k < 4; ++k)
			{
				sim.quatrot[i][j][k] += dt*qdot[k];
			}
			qua_normalize(sim.quatrot[i][j]);
		}
	}
}

void _verlet_rot2(System::simulation& sim)
{
	double dt_half = sim.timestep/2;
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < 3; ++k)
			{
				sim.angmomentum[i][j][k] += dt_half*sim.torque[i][j][k];
			}
		}
	}
}

void _rotate_mol(System::simulation& sim)
{
	#pragma omp target teams distribute parallel for
	for (int i = 0; i < sim.n_types; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			for (int l = 0; l < sim.n_atoms[i]; ++l)
			{
				matrix_mult_vec(sim.rotation_matrix[i][j],sim.position_par_com_init[i][j][l],sim.position_par_com[i][j][l]);
				for (int k = 0; k < 3; ++k)
				{
					sim.position_par_world[i][j][l][k] = sim.position_com[i][j][k] + sim.position_par_com[i][j][l][k];
				}
			}
		}
	}
}

/* Internal Functions end here*/
//
//
//
//
// MAINTAIN DISTANCE
//
//
//
//
//

void integrate_verlet_periodic(System::simulation& sim)
{
	double dt = sim.timestep;
	sim.energy_total = sim.energy_potential;
	
	//First part of velocity verlet
	_verlet_trans1(sim); // Translation P1
	_verlet_rot1(sim); // Rotation P2

	_get_rotation_matrix(sim);
	_rotate_mol(sim);
	//Calling interaction
	sim.energy_potential = 0;
	interact(sim);
	//Done

	//Second part of velocity verlet
	_verlet_trans2(sim); //Translation P2
	_verlet_rot2(sim); // Rotation P2

	sim.time+= dt;
	sim.state++;
}

/*******************************************************************************
 * This function integrates the equation of motion for rigid box conditions
 * This function integrates using the Leapfrog Algorithm 
 *
 * @param sim Simulation being integrated over
 ******************************************************************************/
void integrate_verdet_box(System::simulation& sim)
{
	double dt = sim.timestep;
	sim.energy_total = sim.energy_potential;
	for (int i = 0; i < sim.n_types; ++i)
	{
		sim.energy_kinetic[i] = 0;
		#pragma omp target teams distribute parallel for
		for (int j = 0; j < sim.n_molecules[i]; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < sim.n_dimensions; ++k)
			{
				sim.position_com[i][j][k] += dt*sim.velocity_com[i][j][k] + 0.5*dt*dt*sim.acceleration_com[i][j][k];
				sim.velocity_com[i][j][k] += dt*sim.acceleration_com[i][j][k];
				if(sim.position_com[i][j][k]<0 || sim.position_com[i][j][k] > sim.box_size_limits[k]){
					int num_bounce = (int)(sim.position_com[i][j][k]/sim.box_size_limits[k]);
					double l = std::abs(sim.position_com[i][j][k] - sim.box_size_limits[k]*num_bounce);
					//Getting rid of the signs
					num_bounce *= ((num_bounce < 0)*-1 + (num_bounce>0));
					//Done
					//Implementing reflection
					sim.position_com[i][j][k] = (num_bounce%2 == 0)*l + (num_bounce%2 == 1)*(sim.box_size_limits[k]-l);
					sim.velocity_com[i][j][k] *= ((num_bounce%2==0) + (num_bounce%2 == 1)*(-1));
					//Done
				}
				#pragma omp atomic
				sim.energy_kinetic[i] += sim.velocity_com[i][j][k]*sim.velocity_com[i][j][k];
			}
		}
		sim.energy_kinetic[i] *= sim.mass[i];
		sim.temperature[i] = sim.energy_kinetic[i]/(sim.n_dimensions*sim.n_molecules[i]*BOLTZ_SI);
		sim.energy_total += sim.energy_kinetic[i];
	}
	sim.time+=dt;
	sim.state++;
}