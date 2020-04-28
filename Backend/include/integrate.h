/** @file */ 
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <cmath>
#include <algorithm>
#include "system.h"
#include "constants.h"
#include "interaction.h"
#include "thermostat.h"

/*******************************************************************************
 * This function integrates the equation of motion for periodic boundary conditions
 * This function integrates using the Leapfrog Algorithm 
 *
 * @param sim Simulation being integrated over
 ******************************************************************************/
void integrate_verdet_periodic(System::simulation& sim);

/*******************************************************************************
 * This function integrates the equation of motion for rigid box conditions
 * This function integrates using the Leapfrog Algorithm 
 *
 * @param sim Simulation being integrated over
 ******************************************************************************/
void integrate_verdet_box(System::simulation& sim);
#endif
