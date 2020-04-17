#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <cmath>
#include <algorithm>
#include "system.h"
#include "constants.h"
#include "interaction.h"
#include "thermostat.h"

void integrate_verdet_periodic(System::simulation&);
void integrate_verdet_box(System::simulation&);
#endif
