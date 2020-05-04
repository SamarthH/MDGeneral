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
 * @param type1 Type of molecule 1
 * @param type2 Type of molecule 2
 * @param n1 Molecule index of molecule 1 in position_par_world[type1] array
 * @param n2 Molecule index of molecule 1 in position_par_world[type1] array
 * @param m1 Particle index of particle in molecule 1
 * @param m2 Particle index of particle in molecule 2
 ******************************************************************************/
double distance_periodic(System::simulation& sim, int type1, int n1, int m1, int type2, int n2, int m2);

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
 * Setup of free particle interaction of atoms n1 and n2 of molecules of type1 and type2.
 * This does nothing. The function is empty.
 *
 * @param sim Simulation being used
 * @param type1 First type of particle interacting
 * @param type2 Second type of particle interacting
 * @param k Atom n1 of type1
 * @param l Atom n2 of type2
 ******************************************************************************/
void free_particles(System::simulation& sim,int type1,int type2, int k, int l);

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for periodic boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for periodic boundary conditions between atoms n1 and n2 of molecules of type1 and type2.
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
 * @param k Atom n1 of type1
 * @param l Atom n2 of type2
 ******************************************************************************/
void lj_periodic(System::simulation& sim,int type1,int type2, int k, int l); 

/*******************************************************************************
 * \brief Setup the Lennard-Jones potential for rigid box boundary conditions between two particle types
 * 
 * Setup of Lennard-Jones potential for rigid box boundary conditions between atoms n1 and n2 of molecules of type1 and type2 \n
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
 * @param k Atom n1 of type1
 * @param l Atom n2 of type2
 ******************************************************************************/
void lj_box(System::simulation& sim,int type1,int type2, int k, int l);

#endif