#ifndef AUX_FUNC_H
#define AUX_FUNC_H

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

#include "dynamical_system.h"
#include "auxiliar_functions_vmo.h"

/******************** simple functions ***************************/

// returns x coordinate of selected collinear lagrangian point
double collinear_lagrangian_point_x_coordinate(void *params, int n);

// returns the distance to one of the primaries
double distance_to_primary(double *y, void *params, double t, char *system, int n);

// calculates the index velocity_index of state array y
int integral_to_velocity(double *y, void *params, double t, double constant, char *system, int velocity_index);

// calculates xdot of point in the rorational reference frame
double xdot(double *y, void *params, double constant);

// calculates ydot of point in the rorational reference frame
double ydot(double *y, void *params, double constant);

// returns energy of selected lagrangian point
double lagrangian_point_jacobi_constant(void *params, int n);

// returns the integral of motion numerical value
double motion_integral(double *y, void *params, double t, char *system);

// calculates the jacobi integral of the system
// regularized around the bigger mass 
// must return zero
double motion_integral_regularized_bigger_mass(double *v, void *params);

// calculates the jacobi integral of the system
// regularized around the smaller mass 
// must return zero
double motion_integral_regularized_smaller_mass(double *v, void *params);

// returns the position of the bigger primary
double primary_position_big(void *params);

// returns the position of the smaller primary
double primary_position_small(void *params);

// returns the pseudo potential value 
// for the rotational system
double pseudo_potential(double *y, void *params);

// returns the valur of the derivative on y of
// the pseudo potential for the rotational system
double pseudo_potential_y_derivative(double *y, void *params);

// shows lagrangian points energy on terminal
int show_lagrangian_points_jacobi_constant(void *params);

// switches from regularized variables regarding the bigger
// mass to rotational variables using Levi-Civita transformation
int switch_back_from_regularized_bigger_mass (double *y, double *v, void *params);

// switches from rotational variables to regularized variables
// regarding the bigger mass using Levi-Civita transformation
int switch_to_regularized_bigger_mass (double *v, double *y, void *params);

// switches from regularized variables regarding the smaller
// mass to rotational variables using Levi-Civita transformation
int switch_back_from_regularized_smaller_mass (double *y, double *v, void *params);

// switches from rotational variables to regularized variables
// regarding the smaller mass using Levi-Civita transformation
int switch_to_regularized_smaller_mass (double *v, double *y, void *params);

// defines y as the symmetric of x according to the planar 
// circular and restricted three body problem symmetry 
int symmetry_pcrtbp(double *y, double *x);

// writes lagrangian points positions to file
int write_lagrangian_points_position(void *params);

// writes masses position to file
int write_masses_position(void *params);

/****************** high-level functions *************************/

// calculates an orbit with given initial condition in the rotational system
// using gsl with rk8pd as the integrator of choice
// defines passed arrays with the complete orbit evenly spaced in time
int evolve(double *y, void *params, double tf, double ***orbit, double **time, int *step_number, int *size);

// calculates an orbit with given initial condition in the rotational system
// together with the transition matrix
// using gsl with rk8pd as the integrator of choice
int evolve_extended(double *y, void *params, double tf);

// calculates the poincare map of a given orbit
// and defines it in the passed arrays
int poincare_map(void *params, double ***map, double **energy, double **bissection_times, int *step_number_map, int *size_map, double **orbit, double *time, int step_number);

/******************* zero-velocity curve *************************/

// field for tracing zero velocity curves
// on the rotational frame
int field_zero_velociy_curves_rotational(double t, const double y[], double f[], void *params);

// structure with parameters for root calculation
struct root_params
{
	double mu, J;
};

// function to be solved for root calculation
// with seed close to the third Lagrangian point
double root_function_L3(double x, void *params);

// derivative of function to be solved
// for root calculation
// with seed close to the third Lagrangian point
double root_derivative_L3(double x, void *params);

// both function and derivative to be solved
// for root calculation
// with seed close to the third Lagrangian point
void root_fdf_L3(double x, void *params, double *y, double *dy);

// function to be solved for root calculation
// with seed close to the fourth and fifth Lagrangian points
double root_function_triangular_points(double y, void *params);

// derivative of function to be solved
// for root calculation
// with seed close to the fourth and fifth Lagrangian points
double root_derivative_triangular_points(double y, void *params);

// both function and derivative to be solved
// for root calculation
// with seed close to the fourth and fifth Lagrangian points
void root_fdf_triangular_points(double y, void *params, double *z, double *dz);

// determines the root for zero velocity curve tracing
int root_zvc(double J, void *params, double *y, double x);

/*********************** periodic orbit **************************/

// calculates the initial condition of a periodic orbit
// also updates its motion constant and period values
int find_periodic_orbit_initial_parameters(void *params, double *y0, double tf, double *Jupo, double *period);

// calculates a Lyapunov orbit in the linear region
// of one of the Lagrangian points
int Lyapunov_orbit_initial_condition_linear_region(void *params, double *y0, double delta, double *J, double *period, int n);

// calculates the Lyapunov orbit for given energy
int find_Lyapunov_orbit_parameters(void *params, double Jupo, double *y0upo, double *periodupo, int n);

/************************ eigen system ***************************/

// returns the determinant of a four by four matrix
// given in array form of size sixteen
double determinant(double data[]);

// calculates the eigenvalues and eigenvectors
// of a four by four matrix given in array form
int calculate_eigen_values_and_vectors(double data[], double eig_value_real[], double eig_value_imag[], double eig_vec_real[][4], double eig_vec_imag[][4]);

/******************** special lunar orbits ***********************/

// finds the parameters of the outer stable
// orbit in the lunar realm
int find_outer_lunar_orbit(void *params, double J_orbit, double *y0_orbit, double *period_orbit, char* stability);

// writes position of outer orbits and
// real part of eigenvalues for a range of J
int profile_outer_lunar_orbit(void *params);

// writes position of outer orbits
// without using a continuation method
int profile_outer_lunar_orbit_v2(void *params);

// writes position of inner orbit and
// real part of eigenvalues for a range of J
int profile_inner_lunar_orbit(void *params);

// writes position of inner orbit for selected
// values of J
int profile_inner_lunar_orbit_v2(void *params);

// traces the manifolds of the saddle
// at the outer lunar orbit
int trace_manifolds_outer_saddle(double J_orbit);

// finds the parameters of the unstable periodic
// orbit around the inner lunar orbit
int find_satellite_UPO(double y0[4], double *J, double *period, char* island);

// traces the manifolds of the upo
// around the inner lunar orbit
int trace_manifolds_satellite_UPO(double *J_orbit, char* island);

#endif