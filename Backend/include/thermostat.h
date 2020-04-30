/** @file */ 
#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "system.h"
#include "interaction.h"
#include "constants.h"
#include "integrate.h"
#include "algorithm_constants.h"
#include <cmath>
#include <algorithm>
#include <random>

/*******************************************************************************
 * \brief Initializes the constant arrays of thermostats for speed
 * 
 * The function initializes the arrays for more efficient computation by precomputing the required factors 
 *
 * @param sim Simulation being initialized
 ******************************************************************************/
void initialize_thermostats(System::simulation& sim);

/*******************************************************************************
 * \brief This function calls all the required thermostats
 * 
 * Calls all the thermostats and updates the velocities of particles according to the chosen thermostat.
 *
 * @param sim Simulation being used
 ******************************************************************************/
void call_thermostat(System::simulation& sim);

/*******************************************************************************
 * \brief Call when no thermostat is to be used
 *
 * Does nothing
 *
 * @param sim Simulation being used
 * @param type Type of particle for thermalization
 ******************************************************************************/
void no_thermostat(System::simulation& sim, int type);

/*******************************************************************************
 * \brief Call when anderson thermostat is to be used
 *
 * Applies the anderson thermostat on the required particle types
 * 
 * For this, a particle is given a speed from the gaussian distribution with 0 mean and \f$ \frac{k_{b}*T_{req}}{m} \f$ variance with a probability of \f$ \nu \f$ per unit time.  \n
 *	 \n
 *	Here, the constants vector is as follows: \n
 *	thermostat_const[i][0] = \f$ \nu \f$ \n
 *	thermostat_const[i][1] = \f$ \sqrt{\frac{k_b*T_{req}}{m}} \f$ for the particle type i
 *
 * @param sim Simulation being used
 * @param type Type of particle for thermalization
 ******************************************************************************/
void anderson(System::simulation& sim,int type);

/*******************************************************************************
 * \brief Call when Bussi–Donadio–Parrinello thermostat is to be used
 *
 * Applies the Bussi–Donadio–Parrinello thermostat on the required particle types.
 * Implemented according to the algorithm described in J. Chem. Phys. 126, 014101 (2007); https://doi.org/10.1063/1.2408420 \n
 * On every timestep, we scale the velocities by, \f$ \alpha \f$, where,
 * \f[ \alpha^2 = a + \frac{b}{K}\sum_{i=1}^{N_f} R_i^2 + c \frac{R_1}{\sqrt{K}} \f]
 * Where, \f$ R_i \f$ are gaussian random numbers, K is kinetic energy of system, \f$ N_f \f$ is number of degrees of freedom (translational) \n
 * Also,
 * \f[ a =  e^{-\Delta t/\tau} \f]
 * \f[ b = (1-a)\frac{\bar{K}}{N_f} = (1-a)\sqrt{T_{req}/2}\f] 
 * \f[ c = 2\sqrt{ab} \f]
 * \n
 *	Here, the constants vector is as follows: \n
 *	thermostat_const[i][0] = \f$ tau \f$ (the relaxation time parameter)\n
 *	thermostat_const[i][1] = \f$ a \f$ \n
 *  thermostat_const[i][2] = \f$ b \f$ \n
 *  thermostat_const[i][3] = \f$ c \f$ \n
 *  \n
 *  Here, exponentiation is done by using a fourth order taylor expansion as it is upto 3 times more efficient and more accurate (even for |x| < 10^-4 ) than the cpp implementation
 *
 * @param sim Simulation being used
 * @param type Type of particle for thermalization
 ******************************************************************************/
 void bussi(System::simulation& sim,int type);
#endif