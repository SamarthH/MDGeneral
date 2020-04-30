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
 * @param sim Simulation being used
 * @param type Type of particle for thermalization
 ******************************************************************************/
void anderson(System::simulation& sim,int type);
#endif