#ifndef AUX_GSL_H
#define AUX_GSL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>

#include "dynamical_system.h"

/**
 * this library deals with functions that uses GSL
 * defined variables as arguments
**/

// sets driver for chosen error and control
int set_driver(gsl_odeiv2_driver *d, gsl_odeiv2_system *sys, gsl_odeiv2_step_type *T, char *error, char *control);

// sets the chosen integrator
int set_integrator(const gsl_odeiv2_step_type *T, char *integrator);

// sets the chosen system
int set_system(gsl_odeiv2_system *sys, void *par, char *system);

#endif