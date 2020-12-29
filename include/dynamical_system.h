#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/**
 * ODE fields and (some) jacobians
 * for different forms of the 
 * planar, circular and restricted 
 * three-body problem
**/

// field for rotational frame plus transition matrix
int field_extended(double t, const double y[], double f[], void *params);

// field for Hamiltonian form
int field_hamiltonian(double t, const double y[], double f[], void *params);

// field for inertial frame
int field_inertial(double t, const double y[], double f[], void *params);

// field for regularized system around the bigger mass
// transformed from the rotational frame and
// valid as a linear approximation
// using Levi-Civita transformation
int field_regularized_bigger_mass(double t, const double y[], double f[], void *params);

// field for regularized system around the smaller mass
// transformed from the rotational frame and
// valid as a linear approximation
// using Levi-Civita transformation
int field_regularized_smaller_mass(double t, const double y[], double f[], void *params);

// field for rotational frame
int field_rotational(double t, const double y[], double f[], void *params);

// jacobian for rotational frame
int jacobian_rotational(double t, const double y[], double *dfdy, double dfdt[], void *params);

#endif