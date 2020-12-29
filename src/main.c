#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main_functions.h"

int main(int argc, char **argv)
{
	/********************** Start clock ************************/

	clock_t begin = clock(), end;

	/********************* Mass parameter **********************/

	double mu = 1.215e-2; // Earth-Moon system

	/******************* Declared variables ********************/

	double J;			// Jacobi constant value
	double time;		// integration time
	double ic[4];		// initial condition ic=(x,y,xdot,ydot)
	double x_min;		// minimum x value
	double x_max;		// maximum x value
	double xdot_min;	// minimum xdot value
	double xdot_max;	// maximum xdot value
	int nx;				// number of initial conditions
						// in x range
	int nxdot;			// number of initial conditions
						// in xdot range
	int i; 				// index of Lagrangian point
	int number_orbits;	// number of Lyapunov orbits
	double delta_J;		// difference of the Jacobi constant
						// between the family elements
	int number_points; 	// number of points taken on the
						// Lyapunov orbit
	double time_left;	// integration time for the
						// left branch of the manifolds
	double time_right;  // integration time for the
						// right branch of the manifolds

	/***********************************************************/
	
	////////////////////////////////////////////////////////////
	/*				     	 Useful info					  */
	////////////////////////////////////////////////////////////

	// show_lagrangian_points_jacobi_constant (&mu);

	////////////////////////////////////////////////////////////
	/*	    Location of Lagrangian points and primaries	      */
	////////////////////////////////////////////////////////////

	// write_lagrangian_points_position (&mu);

	// write_masses_position (&mu);

	////////////////////////////////////////////////////////////
	/*					   		Orbit						  */
	////////////////////////////////////////////////////////////

	// J = 3.172;
	// time = 100.0;
	// ic[0] = 1.1;
	// ic[1] = 0.0;
	// ic[2] = 0.0;
	// ic[3] = ydot(ic, &mu, J);	// calculates ydot based on
	// 							// the Jacobi constant value

	// trace_orbit(ic, &mu, time);

	////////////////////////////////////////////////////////////
	/*				     Zero-velocity curve				  */
	////////////////////////////////////////////////////////////

	// J = 3.172;
	// trace_zvc(J, &mu);

	////////////////////////////////////////////////////////////
	/*				  		 Phase Space 					  */
	////////////////////////////////////////////////////////////

	// J = 3.172;
	// time = 2000.0;
	// x_min = -0.5;
	// x_max = 1.11375;
	// xdot_min = -3.0;
	// xdot_max = 3.0;
	// nx = 15;
	// nxdot = 15;

	// draw_phase_space(J, &mu, time, x_min, x_max,
	// 					xdot_min, xdot_max, nx, nxdot); 

	////////////////////////////////////////////////////////////
	/*					  Lyapunov orbit					  */
	////////////////////////////////////////////////////////////

	// J = 3.172;
	// i = 1; 
	// trace_Lyapunov_orbit(&mu, J, i);

	// J = 3.172;
	// i = 2;
	// trace_Lyapunov_orbit(&mu, J, i);

	// i = 1;
	// number_orbits = 10; 
	// delta_J = 0.02;
	// trace_Lyapunov_orbit_family(&mu, number_orbits,
	// 							delta_J, i);

	// i = 2;
	// number_orbits = 9; 
	// delta_J = 0.02;
	// trace_Lyapunov_orbit_family(&mu, number_orbits,
	// 							delta_J, i);
	
	////////////////////////////////////////////////////////////
	/*					 Manifolds tracing					  */
	////////////////////////////////////////////////////////////

	// J = 3.172;
	// i = 1;
	// number_points = 100;
	// time_left = 12.0;
	// time_right = 5.0;
	// trace_manifolds_Lyapunov_orbit(&mu, J,
	// 				time_left, time_right, number_points, i);

	// J = 3.172;
	// i = 2; 
	// number_points = 100;
	// time_left = 8.0;
	// time_right = 8.0;
	// trace_manifolds_Lyapunov_orbit(&mu, J,
	// 				time_left, time_right, number_points, i);

	/********************** Stop clock *************************/

	end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	if (time_spent < 60.0)
	{
		printf("time spent = %1.2e seconds\n", time_spent);
	}
	else if (time_spent < 3600.0)
	{
		printf("time spent = %1.2e minutes\n", time_spent/60.0);
	}
	else
	{
		printf("time spent = %1.2e hours\n", time_spent/3600.0);
	}
	

	/***********************************************************/

	return 0;
}