#ifndef MAIN_FUNC_H
#define MAIN_FUNC_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "auxiliar_functions.h"
#include "auxiliar_functions_gsl.h"
#include "auxiliar_functions_vmo.h"
#include "dynamical_system.h"

// calculates and writes an orbit with given initial 
// condition in the rotational system
// also calculates and writes orbit on poincare map
int trace_orbit(double *y, void *params, double tf);

// calculates and writes the phase space  
// for a given jacobi constant value
int draw_phase_space(double motion_constant, void *params, double tf, double coordinate_min, double coordinate_max, double velocity_min, double velocity_max, int nc, int nv);

// writes the Lyapunov orbit for given energy
int trace_Lyapunov_orbit(void *params, double J, int n);

// writes family of Lyapunov orbits for 
// n-th Lagrangian point
int trace_Lyapunov_orbit_family(void *params, int number_of_orbits, double delta_J, int n);

// traces and writes zero-velocity curve
int trace_zvc(double J, void *params);

// traces and writes the manifolds of a Lyapunov orbit 
int trace_manifolds_Lyapunov_orbit(void *params, double J, double time_left, double time_right, int man_ic_number, int n);

#endif