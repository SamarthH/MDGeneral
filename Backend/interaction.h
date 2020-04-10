#include "system.h"
#include <omp.h>

double distance(simulation*, int a, int b); //Gives the distance between two particles labelled number a,b in the position array

void free_particles_periodic(simulation*);//Free particles with periodic boundary conditions

void free_particles_box(simulation*);//Free particles contained within a box
