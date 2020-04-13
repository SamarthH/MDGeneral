#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <cmath>
#include <algorithm>
#include "system.h"
#include "constants.h"
#include "interaction.h"

void integrate_verdet_periodic(simulation*);
void integrate_verdet_box(simulation*);
#endif
