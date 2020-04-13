#ifndef INTERACTION_H
#define INTERACTION_H

#include "system.h"

double distance(simulation*, int type1, int n1, int type2, int n2); //Gives the distance between two particles labelled number n1,n2 of type1,type2 respectively in the position array

void interact(simulation*); //Calls the required interaction functions

void free_particles(simulation*,int,int);//Free particles

void lj_periodic(simulation*,int,int); //Lennard-Jones potential for periodic boundary conditions and cutoff = half the box size

#endif