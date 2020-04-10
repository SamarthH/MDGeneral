#ifndef INTERACTION_H
#define INTERACTION_H

#include "system.h"

double distance(simulation*, int type1, int n1, int type2, int n2); //Gives the distance between two particles labelled number n1,n2 of type1,type2 respectively in the position array

void free_particles_periodic(simulation*);//Free particles with periodic boundary conditions

void free_particles_box(simulation*);//Free particles contained within a box

#endif