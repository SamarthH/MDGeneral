#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "system.h"
#include "interaction.h"
#include "constants.h"
#include "integrate.h"
#include "algorithm_constants.h"
#include <cmath>
#include <algorithm>
#include <random>

void call_thermostat(simulation*); //Call this to get thermostat effects

void no_thermostat(simulation*, int);

void anderson(simulation*,int);
#endif