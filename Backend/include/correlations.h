/** @file */
#ifndef CORRELATIONS_H
#define CORRELATIONS_H

#include "system.h"

/*******************************************************
 * \brief This function initializes the correlation arrays by inputing intial values
 *
 * This function fills in the initial values (t=0) of the vectors over which correlation is to be found. \n
 * This should be called only after initializing the initial vectors over which correlation is to be found.
 *
 * @param sim Simulation being used
*/
void initialize_correlations(System::simulation& sim);

/*******************************************************
 * \brief Fills in the correlation vectors for this timestep
 *
 * @param sim Simulation being used
*/
void correlate(System::simulation& sim);


#endif