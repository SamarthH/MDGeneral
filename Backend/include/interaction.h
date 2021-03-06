/** @file */ 
#ifndef INTERACTION_H
#define INTERACTION_H

#include "system.h"

/*******************************************************************************
 * \brief Initializes the constant arrays for interactions for speed
 * 
 * The function initializes the arrays for more efficient computation by precomputing the required factors 
 *
 * @param sim Simulation being initialized
 ******************************************************************************/
void initialize_interactions(System::simulation& sim);

/*******************************************************************************
 * \brief Returns the distance between two particles for periodic boundary conditions
 * 
 * Returns the distance between two particles labelled number n1,n2 of type1,type2 
 * respectively in the position array assuming periodic boundary conditions
 *
 * @param sim Simulation being used
 * @param type1 Particle type of particle 1
 * @param type2 Particle type of particle 2
 * @param n1 Particle index of particle 1 in position[type1] array
 * @param n2 Particle index of particle 2 in position[type2] array
 ******************************************************************************/
double distance_periodic(System::simulation& sim, int type1, int n1, int type2, int n2);

/*******************************************************************************
 * \brief This function calls all the required interaction functions between the particles
 * 
 * Calls all the interactiosn and updates acceleration and energy arrays
 *
 * @param sim Simulation being used
 ******************************************************************************/
void interact(System::simulation& sim);

/*******************************************************************************
 * \brief Setup the free particle interaction between two particle types
 * 
 * Setup of free particle interaction between particle types type1 and type2.
 * This does nothing. The function is empty.
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void free_particles(System::simulation& sim,int type1,int type2);

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for periodic boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for periodic boundary conditions between particle types type1 and type2.
 * This requires the cutoff <= half the box size because it only checks the nearest images. \n
 * \n
 * Potential : \f[ U(r) = 4\epsilon*[(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6] - U(r_cut) \f] \n
 * Force : \f[ F(r) = 24\epsilon*\frac{\sigma^6}{r^7}*[2(\frac{\sigma}{r})^6 - 1] \vu{r} \f] \n
 * Here, the constants vector is as follows: \n
 * interaction_const[i][j][0] = \f$ \epsilon \f$ \n
 * interaction_const[i][j][1] = \f$ \sigma \f$ \n
 * interaction_const[i][j][2] = Cutoff radius/distance (r_cut) \n
 * interaction_const[i][j][3] = Truncated Potential (etrunc) \n
 * interaction_const[i][j][4] = \f$ \sigma^6 \f$ \n
 * interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius) \n
 * \n
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void lj_periodic(System::simulation& sim,int type1,int type2); 

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for rigid box boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for rigid box boundary conditions between particle types type1 and type2. \n
 * \n
 * Potential : \f[ U(r) = 4\epsilon*[(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6] - U(r_cut) \f] \n
 * Force : \f[ F(r) = 24\epsilon*\frac{\sigma^6}{r^7}*[2(\frac{\sigma}{r})^6 - 1] \vu{r} \f] \n
 * Here, the constants vector is as follows: \n
 * interaction_const[i][j][0] = \f$ \epsilon \f$ \n
 * interaction_const[i][j][1] = \f$ \sigma \f$ \n
 * interaction_const[i][j][2] = Cutoff radius/distance (r_cut) \n
 * interaction_const[i][j][3] = Truncated Potential (etrunc) \n
 * interaction_const[i][j][4] = \f$ \sigma^6 \f$ \n
 * interaction_const[i][j][5] = Tail Energy (assuming constant distribution outside cutoff radius) \n
 * \n
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 ******************************************************************************/
void lj_box(System::simulation& sim,int type1,int type2);

#endif