#ifndef INTERACTION_H
#define INTERACTION_H

#include "system.h"

/*Initializes the constant arrays for interactions for speed*/
void initialize_interactions(System::simulation& sim);

/*Gives the distance between two particles labelled number n1,n2 of type1,type2 respectively in the position array*/
double distance(System::simulation&, int type1, int n1, int type2, int n2);

/*Calls the required interaction functions*/
void interact(System::simulation&);

/*Free particles*/
void free_particles(System::simulation&,int,int);

/*Lennard-Jones potential for periodic boundary conditions and cutoff <= half the box size*/
void lj_periodic(System::simulation&,int,int); 

/*Lennard-Jones potential for periodic boundary conditions and cutoff < box size*/
void lj_box(System::simulation&,int,int);

#endif