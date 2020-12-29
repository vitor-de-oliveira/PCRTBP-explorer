#include "auxiliar_functions.h"

double collinear_lagrangian_point_x_coordinate(void *params, int n)
{
	double result;
	double mu = *(double *)params;
	double mu2, mu3, mu4;
	double x2, x3, x4, x5;
	double a0, a1, a2, a3, a4, a5;
	double x, x_new, p, p_derivative, delta;
	mu2 = mu * mu;
	mu3 = mu * mu * mu;
	mu4 = mu * mu * mu * mu;
	if (n == 1)
	{
		a0 = -1. + 3. * mu - 3. * mu2 + 2. * mu3;
		a1 = 2. - 4. * mu + 5. * mu2 - 2. * mu3 + mu4;
		a2 = -1. + 4. * mu - 6. * mu2 + 4. * mu3;
		a3 = 1. - 6. * mu + 6. * mu2;
		a4 = -2. + 4. * mu;
		a5 = 1.;
	}
	else if (n == 2)
	{
		a0 = -1. + 3. * mu - 3. * mu2;
		a1 = 2. - 4. * mu + mu2 - 2. * mu3 + mu4;
		a2 = -1. + 2. * mu - 6. * mu2 + 4. * mu3;
		a3 = 1. - 6. * mu + 6. * mu2;
		a4 = -2. + 4. * mu;
		a5 = 1.;
	}
	x = (1. - mu) - pow((mu / 3.), (1. / 3.));
	do
	{
		x2 = x * x;
		x3 = x * x * x;
		x4 = x * x * x * x;
		x5 = x * x * x * x * x;
		p = a0 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4 + a5 * x5;
		p_derivative = a1 + 2. * a2 * x + 3. * a3 * x2 + 4. * a4 * x3 + 5. * a5 * x4;
		x_new = x - p / p_derivative;
		delta = fabs(x_new - x);
		x = x_new;
	} while (delta > 1e-15);
	result = x_new;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double frac_mu = mu_2 / mu_1;
	double beta = -(7. / 12.) * frac_mu + (7. / 12.) * frac_mu * frac_mu - (13223. / 20736.) * frac_mu * frac_mu * frac_mu;
	if (n == 3) result = mu_1 - (2. + beta);
	if (n != 1 && n != 2 && n != 3)
	{
		printf("Warning: invalid integer for\ncollinear lagrangian function\n");
		exit(2);
	}
	return result;
}

double distance_to_primary(double *y, void *params, double t, char *system, int n)
{
	double result, mu, mu_1, mu_2;
	mu = *(double *)params;
	mu_1 = 1.0 - mu;
	mu_2 = mu;
	if (n != 1 && n != 2)
	{
		printf("Warning: invalid primary id\nfor distance calculation\n");
		exit(2);
	}
	if (strcmp(system, "rotational") == 0 || 
		strcmp(system, "extended") == 0 || 
		strcmp(system, "regularized") == 0 ||
		strcmp(system, "hamiltonian") == 0)
	{
		if (n == 1)
		{
			result = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
		}
		else
		{
			result = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
		}
	}
	else if (strcmp(system, "inertial") == 0)
	{
		if (n == 1)
		{
			result = sqrt((y[0] + mu_2 * cos(t)) * (y[0] + mu_2 * cos(t)) + (y[1] + mu_2 * sin(t)) * (y[1] + mu_2 * sin(t)));
		}
		else
		{
			result = sqrt((y[0] - mu_1 * cos(t)) * (y[0] - mu_1 * cos(t)) + (y[1] - mu_1 * sin(t)) * (y[1] - mu_1 * sin(t)));
		}
	}
	else
	{
		printf("Warning: unknown system type\n");
		exit(2);
	}
	return result;
}

int integral_to_velocity(double *y, void *params, double t, double constant, char *system, int velocity_index)
{
	double mu = *(double *)params;
	double mu_1 = 1.0 - mu;
	double mu_2 = mu;
	double r_1 = distance_to_primary (y, params, t, system, 1);
	double r_2 = distance_to_primary (y, params, t, system, 2);
	double omega;
	if (velocity_index != 2 && velocity_index != 3)
	{
		printf("Warning: invalid velocity index for\nintegral to variable calculation\n");
		exit(2);
	}
	if (strcmp(system, "rotational") == 0 || 
		strcmp(system, "extended") == 0 || 
		strcmp(system, "regularized") == 0)
	{
		omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
		if (velocity_index == 2)
		{
			y[2] = sqrt(2.0 * omega - constant - y[3] * y[3]);
		}
		else if (velocity_index == 3)
		{
			y[3] = sqrt(2.0 * omega - constant - y[2] * y[2]);
		}
	}
	else if (strcmp(system, "hamiltonian") == 0)
	{
		omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
		if (velocity_index == 2)
		{
			y[2] = sqrt(2. * (omega + constant) - (y[3] - y[0]) * (y[3] - y[0]) ) - y[1];
		}
		else if (velocity_index == 3)
		{
			y[3] = sqrt(2. * (omega + constant) - (y[2] + y[1]) * (y[2] + y[1]) ) + y[0];
		}
	}
	else if (strcmp(system, "inertial") == 0)
	{
		omega = mu_1 / r_1 + mu_2 / r_2;
		if (velocity_index == 2)
		{
			y[2] = -1.0 * y[1] + sqrt(y[1] * y[1] - y[3] * y[3] + 2.0 * y[0] * y[3] + 2.0 * omega - constant);
		}
		else if (velocity_index == 3)
		{
			y[3] = y[0] + sqrt(y[0] * y[0] - y[2] * y[2] - 2.0 * y[1] * y[2] + 2.0 * omega - constant);
		}
	}
	else
	{
		printf("Warning: unknown system type for\nintegral to velocity calculation\n");
		exit(2);
	}
	return 0;
}

double xdot(double *y, void *params, double constant)
{
	double mu = *(double *)params;
	double mu_1 = 1.0 - mu;
	double mu_2 = mu;
	double r_1 = distance_to_primary (y, params, 0.0, "rotational", 1);
	double r_2 = distance_to_primary (y, params, 0.0, "rotational", 2);
	double omega;
	omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
	return sqrt(2.0 * omega - constant - y[3] * y[3]);
}

double ydot(double *y, void *params, double constant)
{
	double mu = *(double *)params;
	double mu_1 = 1.0 - mu;
	double mu_2 = mu;
	double r_1 = distance_to_primary (y, params, 0.0, "rotational", 1);
	double r_2 = distance_to_primary (y, params, 0.0, "rotational", 2);
	double omega;
	omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
	return sqrt(2.0 * omega - constant - y[2] * y[2]);
}

double lagrangian_point_jacobi_constant(void *params, int n)
{
	double result, mu, y[4];
	mu = *(double *)params;
	if (n == 1)
	{
		y[0] = collinear_lagrangian_point_x_coordinate(&mu, 1);
		y[1] = 0.0;
		y[2] = 0.0;
		y[3] = 0.0;
	}
	else if (n == 2)
	{
		y[0] = collinear_lagrangian_point_x_coordinate(&mu, 2);
		y[1] = 0.0;
		y[2] = 0.0;
		y[3] = 0.0;
	}
	else if (n == 3)
	{
		y[0] = collinear_lagrangian_point_x_coordinate(&mu, 3);
		y[1] = 0.0;
		y[2] = 0.0;
		y[3] = 0.0;
	}
	else if (n == 4)
	{
		y[0] = 0.5 - mu;
		y[1] = sqrt(3.0) / 2.0;
		y[2] = 0.0;
		y[3] = 0.0;
	}
	else if (n == 5)
	{
		y[0] = 0.5 - mu;
		y[1] = -sqrt(3.0) / 2.0;
		y[2] = 0.0;
		y[3] = 0.0;
	}
	else
	{
		printf("Warning: invalid integer for\nlagrangian energy calculation\n");
		exit(2);
	}
	result = motion_integral(y, &mu, 0.0, "rotational");
	return result;
}

double motion_integral(double *y, void *params, double t, char *system)
{
	double mu = *(double *)params;
	double mu_1 = 1.0 - mu;
	double mu_2 = mu;
	double r_1 = distance_to_primary (y, params, t, system, 1);
	double r_2 = distance_to_primary (y, params, t, system, 2);
	double omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
	double v_squared = y[2] * y[2] + y[3] * y[3];
	double arg_1 = 0.5 * (y[2] + y[1]) * (y[2] + y[1]);
	double arg_2 = 0.5 * (y[3] - y[0]) * (y[3] - y[0]);
	double omega_inertial = mu_1 / r_1 + mu_2 / r_2;
	double arg_1_inertial = y[2] * y[2] + y[3] * y[3];
	double arg_2_inertial = y[0] * y[3] - y[1] * y[2];
	if (strcmp(system, "rotational") == 0 || 
		strcmp(system, "extended") == 0 || 
		strcmp(system, "regularized") == 0)
	{
		return 2.0 * omega - v_squared;
	}
	else if (strcmp(system, "hamiltonian") == 0)
	{
		return arg_1 + arg_2 - omega;
	}
	else if (strcmp(system, "inertial") == 0)
	{
		return -1.0 * arg_1_inertial + 2.0 * arg_2_inertial + 2.0 * omega_inertial;
	}
	else
	{
		printf("Warning: unknown system type\nfor motion integral calculation\n");
		exit(2);
	}
}

double motion_integral_regularized_bigger_mass(double *v, void *params)
{
	double *par = (double *)params;
	double mu = par[0];
	double J = par[1];
	double C = J + mu * (1.0 - mu);
	double u_squared = v[0] * v[0];
	double v_squared = v[1] * v[1];
	double arg_1 = (u_squared + v_squared) * (u_squared + v_squared);
	double arg_2 = - 2.0 * mu * (v[0] + v[1]) * (v[0] - v[1]);
	double arg_3 = mu - C;
	double arg_4 = 2.0 * mu / sqrt(1.0 + (u_squared + v_squared) * (u_squared + v_squared) - 2.0 * (v[0] + v[1]) * (v[0] - v[1]));
	double V = 4.0 * (1.0 - mu) + 2.0 * (v[0] * v[0] + v[1] * v[1]) * (arg_1 + arg_2 + arg_3 + arg_4);
	return v[2] * v[2] + v[3] * v[3] - 2.0 * V;
}

double motion_integral_regularized_smaller_mass(double *v, void *params)
{
	double *par = (double *)params;
	double mu = par[0];
	double J = par[1];
	double C = J + mu * (1.0 - mu);
	double u_squared = v[0] * v[0];
	double v_squared = v[1] * v[1];
	double arg_1 = (u_squared + v_squared) * (u_squared + v_squared);
	double arg_2 = 2.0 * (1.0 - mu) * (v[0] + v[1]) * (v[0] - v[1]);
	double arg_3 = 1.0 - mu - C;
	double arg_4 = 2.0 * (1.0 - mu) / sqrt(1.0 + (u_squared + v_squared) * (u_squared + v_squared) + 2.0 * (v[0] + v[1]) * (v[0] - v[1]));
	double V = 4.0 * mu + 2.0 * (v[0] * v[0] + v[1] * v[1]) * (arg_1 + arg_2 + arg_3 + arg_4);
	return v[2] * v[2] + v[3] * v[3] - 2.0 * V;
}

double primary_position_big(void *params)
{
	double mu = *(double *)params;
	return -mu;
}

double primary_position_small(void *params)
{
	double mu = *(double *)params;
	return 1.0 - mu;
}

double pseudo_potential(double *y, void *params)
{
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = distance_to_primary(y, &mu, 0.0, "rotational", 1);
	double r_2 = distance_to_primary(y, &mu, 0.0, "rotational", 2);
	double omega = 0.5 * (y[0] * y[0] + y[1] * y[1]) + mu_1 / r_1 + mu_2 / r_2;
	return omega;
}

double pseudo_potential_y_derivative(double *y, void *params)
{
	double mu = *(double *)params;
	double r_1 = distance_to_primary(y, &mu, 0.0, "rotational", 1);
	double r_2 = distance_to_primary(y, &mu, 0.0, "rotational", 2);
	double omega_dy = y[1] * (1. - (1. - mu) / (r_1*r_1*r_1) - mu / (r_2*r_2*r_2));
	return omega_dy;
}

int show_lagrangian_points_jacobi_constant(void *params)
{
	double mu, J_L1, J_L2, J_L3, J_L4, J_L5;
	mu = *(double *)params;
	J_L1 = lagrangian_point_jacobi_constant(&mu, 1);
	J_L2 = lagrangian_point_jacobi_constant(&mu, 2);
	J_L3 = lagrangian_point_jacobi_constant(&mu, 3);
	J_L4 = lagrangian_point_jacobi_constant(&mu, 4);
	J_L5 = lagrangian_point_jacobi_constant(&mu, 5);
	printf("J1: %1.12e\nJ2: %1.12e\nJ3: %1.12e\nJ4: %1.12e\nJ5: %1.12e\n", J_L1, J_L2, J_L3, J_L4, J_L5);
	return 0;
}

int switch_back_from_regularized_bigger_mass (double *y, double *v, void *params)
{
	double mu = *(double *)params;
	y[0] = (v[0] - v[1]) * (v[0] + v[1]) - mu;
	y[1] = 2.0 * v[0] * v[1];
	y[2] = 0.5 * (v[0] * v[2] - v[1] * v[3]) / (v[0] * v[0] + v[1] * v[1]);
	y[3] = 0.5 * (v[1] * v[2] + v[0] * v[3]) / (v[0] * v[0] + v[1] * v[1]);
	return 0;
}

int switch_to_regularized_bigger_mass (double *v, double *y, void *params)
{
	double mu = *(double *)params;
	double b, c, delta, u_squared;
	b = -1.0 * (mu + y[0]);
	delta = b * b + y[1] * y[1];
	c = - 0.25 * y[1] * y[1];
	if (b < 0)
	{
		u_squared = -0.5 * b + 0.5 * sqrt(delta);
	}
	else
	{
		u_squared = 2.0 * c / (-1.0 * b - sqrt(delta));
	}
	v[0] = sqrt(u_squared);
	v[1] = y[1] / (2.0 * v[0]);
	v[2] = 2.0 * v[0] * y[2] + 2.0 * v[1] * y[3];
	v[3] = 2.0 * v[0] * y[3] - 2.0 * v[1] * y[2];
	return 0;
}

int switch_back_from_regularized_smaller_mass (double *y, double *v, void *params)
{
	double mu = *(double *)params;
	y[0] = 1.0 - mu + (v[0] - v[1]) * (v[0] + v[1]);
	y[1] = 2.0 * v[0] * v[1];
	y[2] = 0.5 * (v[0] * v[2] - v[1] * v[3]) / (v[0] * v[0] + v[1] * v[1]);
	y[3] = 0.5 * (v[1] * v[2] + v[0] * v[3]) / (v[0] * v[0] + v[1] * v[1]);
	
	// double mu_double = *(double *)params;
	// long double mu = (long double) mu_double;
	// // printf("mu = %Le\n", mu);
	// long double y_local[4], v_local[4];
	// for (int i = 0; i < 4; i++)
	// {
	// 	v_local[i] = (long double) v[i];
	// }
	// y_local[0] = 1.0 - mu + (v_local[0] - v_local[1]) * (v_local[0] + v_local[1]);
	// y_local[1] = 2.0 * v_local[0] * v_local[1];
	// y_local[2] = 0.5 * (v_local[0] * v_local[2] - v_local[1] * v_local[3]) / (v_local[0] * v_local[0] + v_local[1] * v_local[1]);
	// y_local[3] = 0.5 * (v_local[1] * v_local[2] + v_local[0] * v_local[3]) / (v_local[0] * v_local[0] + v_local[1] * v_local[1]);
	// for (int i = 0; i < 4; i++)
	// {
	// 	y[i] = (double) y_local[i];
	// }
	return 0;
}

int switch_to_regularized_smaller_mass (double *v, double *y, void *params)
{
	double mu = *(double *)params;
	double b, c, delta, u_squared;
	b = 1.0 - mu - y[0];
	delta = b * b + y[1] * y[1];
	c = -0.25 * y[1] * y[1]; 
	if (b < 0)
	{
		u_squared = -0.5 * b + 0.5 * sqrt(delta);
	}
	else
	{
		u_squared = 2.0 * c / (-1.0 * b - sqrt(delta));
	}
	v[0] = sqrt(u_squared);
	v[1] = y[1] / (2.0 * v[0]);
	v[2] = 2.0 * v[0] * y[2] + 2.0 * v[1] * y[3];
	v[3] = 2.0 * v[0] * y[3] - 2.0 * v[1] * y[2];

	// double mu_double = *(double *)params;
	// long double mu = (long double) mu_double;
	// // printf("mu = %Le\n", mu);
	// long double b, c, delta, u_squared;
	// long double y_local[4], v_local[4];
	// for (int i = 0; i < 4; i++)
	// {
	// 	y_local[i] = (long double) y[i];
	// }
	// // printf("%Le\n", y_local[0]);
	// // printf("%le\n", y[0]);
	// b = 1.0 - mu - y_local[0];
	// delta = b * b + y_local[1] * y_local[1];
	// c = -0.25 * y_local[1] * y_local[1]; 
	// if (b < 0)
	// {
	// 	u_squared = -0.5 * b + 0.5 * sqrtl(delta);
	// }
	// else
	// {
	// 	u_squared = 2.0 * c / (-1.0 * b - sqrtl(delta));
	// }
	// // printf("mu = %Le\nb = %Le\ndelta = %Le\nc = %Le\nu_squared = %Le\n", mu, b, delta, c, u_squared);
	// v_local[0] = sqrtl(u_squared);
	// v_local[1] = y_local[1] / (2.0 * v_local[0]);
	// v_local[2] = 2.0 * v_local[0] * y_local[2] + 2.0 * v_local[1] * y_local[3];
	// v_local[3] = 2.0 * v_local[0] * y_local[3] - 2.0 * v_local[1] * y_local[2];
	// for (int i = 0; i < 4; i++)
	// {
	// 	v[i] = (double) v_local[i];
	// 	// printf("%Le\n", v_local[i]);
	// 	// printf("%le\n", v[i]);
	// }
	return 0;
}

int symmetry_pcrtbp(double *y, double *x)
{
	y[0] = x[0];
	y[1] = -x[1];
	y[2] = -x[2];
	y[3] = x[3];
	return 0;
}

int write_lagrangian_points_position(void *params)
{
	FILE *out = fopen("results/location/lagrangian_points_position.dat", "w");
	double mu, lagrange[5][2];
	mu = *(double *)params;
	lagrange[0][0] = collinear_lagrangian_point_x_coordinate(&mu, 1);
	lagrange[0][1] = 0.0;
	lagrange[1][0] = collinear_lagrangian_point_x_coordinate(&mu, 2);
	lagrange[1][1] = 0.0;
	lagrange[2][0] = collinear_lagrangian_point_x_coordinate(&mu, 3);
	lagrange[2][1] = 0.0;
	lagrange[3][0] = 0.5 - mu;
	lagrange[3][1] = sqrt(3.0) / 2.0;
	lagrange[4][0] = 0.5 - mu;
	lagrange[4][1] = -sqrt(3.0) / 2.0;
	for (int i = 0; i < 5; i++)
	{
		fprintf(out, "%1.16e %1.16e\n", lagrange[i][0], lagrange[i][1]);
	}
	fclose(out);
	printf("Data written in folder results/location\n");
	return 0;
}

int write_masses_position(void *params)
{
	FILE *mas = fopen("results/location/masses_position.dat", "w");
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	fprintf(mas, "%1.16e %1.16e\n", -mu_2, 0.0);
	fprintf(mas, "%1.16e %1.16e\n", mu_1, 0.0);
	fclose(mas);
	printf("Data written in folder results/location\n");
	return 0;
}

int evolve(double *y, void *params, double tf, double ***orbit, double **time, int *step_number, int *size)
{
	// declare variables
	int step_counter, space_orbit, status, system_dimension;
	double mu = *(double *)params, t, t0, h, time_step;
	double const_check, box, constant_error_tolerance;
	double m1_distance, regularization_region_big, collision_radius_big;
	double m2_distance, regularization_region_small, collision_radius_small;
	double const_init = motion_integral(y, &mu, 0.0, "rotational");
	double par_reg[2] = {mu, const_init};
	double v[4], dummy_v[4], dummy_t;
	bool collision_check = false;
	bool inside_regular_big = false;
	bool inside_regular_small = false;

	// initialiaze control variables
	box = 1e3;
	collision_radius_big = 1.66e-2;
	collision_radius_small = 4.52e-3;
	constant_error_tolerance = 1e-8;
	h = 1e-3 * sign(tf);
	regularization_region_big = 3.67e-2;
	regularization_region_small = 1e-2;
	time_step = 1e-3 * sign(tf);
	system_dimension = 4;
	t0 = 0.0;

	// set principal driver
	gsl_odeiv2_system sys = 
	{field_rotational, jacobian_rotational, system_dimension, &mu};
	gsl_odeiv2_driver *d = 
	gsl_odeiv2_driver_alloc_standard_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-15, 0.0, 0.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d, 1e-3);
	gsl_odeiv2_driver_set_hmin(d, 1e-11);

	// set driver for regularized region around the bigger primary
	gsl_odeiv2_system sys_reg_big =
	{field_regularized_bigger_mass, NULL, system_dimension, par_reg};
	gsl_odeiv2_driver *d_reg_big =
	gsl_odeiv2_driver_alloc_standard_new(&sys_reg_big, gsl_odeiv2_step_rk8pd, h, 1e-13, 0.0, 0.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d_reg_big, 1e-3);
	gsl_odeiv2_driver_set_hmin(d_reg_big, 1e-11);

	// set driver for regularized region around the smaller primary
	gsl_odeiv2_system sys_reg_small =
	{field_regularized_smaller_mass, NULL, system_dimension, par_reg};
	gsl_odeiv2_driver *d_reg_small =
	gsl_odeiv2_driver_alloc_standard_new(&sys_reg_small, gsl_odeiv2_step_rk8pd, h, 1e-13, 0.0, 0.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d_reg_small, 1e-3);
	gsl_odeiv2_driver_set_hmin(d_reg_small, 1e-11);
	
	// allocate memory and initializes exit data
	t = t0;
	space_orbit = 1000;
	if (time != NULL)
	{
		alloc_1d_double(time, space_orbit);
		(*time)[0] = t;
	}
	if (orbit != NULL)
	{
		alloc_2d_double(orbit, space_orbit, system_dimension);
		copy((*orbit)[0], y, system_dimension);
	}

	// initialize counter
	step_counter = 1;

	// FILE *tst = fopen("results/orbit/regularized_energy.dat","w");

	// start orbit evolution
	while (fabs(t) < fabs(tf))
	{
		// measures distance to primaries
		m1_distance = distance_to_primary(y, &mu, t, "rotational", 1);
		m2_distance = distance_to_primary(y, &mu, t, "rotational", 2);

		// calculates time evolution
		if (m1_distance > regularization_region_big &&
			m2_distance > regularization_region_small)
		{
			inside_regular_big = false;
			inside_regular_small = false;
			if (fabs(tf - t) < fabs(time_step))
			{
				time_step = tf - t;
			}
			gsl_odeiv2_driver_reset_hstart(d_reg_big, h);
			gsl_odeiv2_driver_reset_hstart(d_reg_small, h);
			status = 
			gsl_odeiv2_driver_apply (d, &t, t + time_step, y);
		}
		else if (m1_distance <= regularization_region_big)
		{
			gsl_odeiv2_driver_reset_hstart(d, h);
			gsl_odeiv2_driver_reset_hstart(d_reg_small, h);
			inside_regular_small = false;
			if (inside_regular_big == false)
			{
				inside_regular_big = true;
				switch_to_regularized_bigger_mass (v, y, &mu);
			}
			if (fabs(tf - t) < fabs(time_step))
			{
				time_step = (tf - t) / (4.0 * (v[0]*v[0] + v[1]*v[1]));
			}
			dummy_t = t;
			copy(dummy_v,v,4);
			status = 
			gsl_odeiv2_driver_apply (d_reg_big, &t, t + time_step, v);
			switch_back_from_regularized_bigger_mass (y, v, &mu);
			t = dummy_t + 2.0 * time_step * (v[0]*v[0] + v[1]*v[1] + dummy_v[0]*dummy_v[0] + dummy_v[1]*dummy_v[1]);
		}
		else
		{
			gsl_odeiv2_driver_reset_hstart(d, h);
			gsl_odeiv2_driver_reset_hstart(d_reg_big, h);
			inside_regular_big = false;
			if (inside_regular_small == false)
			{
				inside_regular_small = true;
				switch_to_regularized_smaller_mass (v, y, &mu);
			}

			/* LOOK FURTHER INTO THIS 11.AUG.2020 */

			if (fabs(tf - t) < fabs(time_step))
			{
				time_step = (tf - t) / (4.0 * (v[0]*v[0] + v[1]*v[1]));
			}
			// if (fabs(tf - t) < fabs(time_step))
			// {
			// 	time_step = 1e-4;
			// }
			// if (fabs(tf - t) < 1e-5)
			// {
			// 	time_step = 1e-8;
			// }

			dummy_t = t;
			copy(dummy_v,v,4);
			status = 
			gsl_odeiv2_driver_apply (d_reg_small, &t, t + time_step, v);
			switch_back_from_regularized_smaller_mass (y, v, &mu);
			t = dummy_t + 2.0 * time_step * (v[0]*v[0] + v[1]*v[1] + dummy_v[0]*dummy_v[0] + dummy_v[1]*dummy_v[1]);
		}
		
		// check if integration was successfull
		if (status != GSL_SUCCESS)
		{
			printf("Warning: %s\n", gsl_strerror(status));
			step_counter--;
			break;
		}

		// check for collision with primaries
		if (collision_check == true)
		{
			if (m1_distance < collision_radius_big)
			{
				// printf("Warning: orbit collided with bigger mass\n");
				step_counter--;
				break;	
			}
			else if (m2_distance < collision_radius_small)
			{
				// printf("Warning: orbit collided with smaller mass\n");
				step_counter--;
				break;	
			}
		}

		// check if orbit reaches box limit
		if (fabs(y[0]) > box || fabs(y[1]) > box)
		{
			printf("Warning: box limit reached\n");
			step_counter--;
			break;
		}

		// check if error in motion constant is acceptable
		if (inside_regular_big == true)
		{
			const_check = motion_integral_regularized_bigger_mass (v, par_reg);
			if (fabs(const_check) > constant_error_tolerance)
			{
				printf("Warning: minimun constant precision reached at regularized system around the bigger mass\n");
				step_counter--;
				break;
			}
		}
		else if (inside_regular_small == true)
		{
			const_check = motion_integral_regularized_smaller_mass (v, par_reg);
			// fprintf(tst, "%1.15e %1.15e\n", t, fabs(const_check));
			if (fabs(const_check) > constant_error_tolerance)
			{
				printf("Warning: minimun constant precision reached at regularized system around the smaller mass\n");
				step_counter--;
				break;
			}
		}
		else
		{
			const_check = motion_integral(y, &mu, t, "rotational");
			if (fabs(const_init - const_check) > constant_error_tolerance)
			{
				printf("Warning: minimum constant precision reached\n");
				step_counter--;
				break;
			}
		}
		
		// allocate more memory space to exit data if necessary
		if (step_counter == space_orbit)
		{
			space_orbit += 1000;
			if (time != NULL)
			{
				realloc_1d_double(time, space_orbit);
			}
			if (orbit != NULL)
			{
				realloc_2d_double(orbit, space_orbit, system_dimension, step_counter);
			}
		}

		// update exit data
		if (time != NULL)
		{
			(*time)[step_counter] = t;
		}
		if (orbit != NULL)
		{
			copy((*orbit)[step_counter], y, system_dimension);
		}

		// update counter
		step_counter++;
	}

	// define exit data size
	if (step_number != NULL)
	{
		*step_number = step_counter;
	}
	if (size != NULL)
	{
		*size = space_orbit;
	}

	// free memory
	gsl_odeiv2_driver_free(d);
	gsl_odeiv2_driver_free(d_reg_big);
	gsl_odeiv2_driver_free(d_reg_small);

	// fclose(tst);

	// printf("t = %1.10f\n", t);

	// printf("time_step = %1.5f\n", time_step);

	return 0;
}

int evolve_extended(double *y, void *params, double tf)
{
	// declare variables
	int status, system_dimension;
	double mu = *(double *)params, t, t0, h, time_step;
	double const_check, constant_error_tolerance;
	double const_init = motion_integral(y, &mu, 0.0, "extended");

	// initialiaze control variables
	constant_error_tolerance = 1e-8;
	h = 1e-3 * sign(tf);
	time_step = 1e-3 * sign(tf);
	system_dimension = 20;
	t0 = 0.0;

	// set driver
	gsl_odeiv2_system sys = 
	{field_extended, NULL, system_dimension, &mu};
	gsl_odeiv2_driver *d = 
	gsl_odeiv2_driver_alloc_standard_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-13, 0.0, 0.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d, 1e-3);
	gsl_odeiv2_driver_set_hmin(d, 1e-11);

	// start orbit evolution
	t = t0;
	while (fabs(t) < fabs(tf))
	{
		// guarantee the calculation will end at final time
		if (fabs(tf - t) < fabs(time_step))
		{
			time_step = tf - t;
		}

		// calculates time evolution
		status = 
		gsl_odeiv2_driver_apply (d, &t, t + time_step, y);

		// check if integration was successfull
		if (status != GSL_SUCCESS)
		{
			printf("Warning: %s\n", gsl_strerror(status));
			break;
		}

		// check if error in motion constant is acceptable
		const_check = motion_integral(y, &mu, t, "extended");
		if (fabs(const_init - const_check) > constant_error_tolerance)
		{
			printf("Warning: minimum constant precision reached\n");
			break;
		}
	}

	// free memory
	gsl_odeiv2_driver_free(d);

	return 0;
}

int poincare_map(void *params, double ***map, double **energy, double **bissection_times, int *step_number_map, int *size_map, double **orbit, double *time, int step_number)
{
	// declare variables
	int space_map, system_dimension;
	int time_counter, bissection_loop_control, bissection_loop_max;
	int map_coordinate, map_velocity;
	int cross_coordinate, cross_velocity;
	double y_bis[4], y_bis_temp[4];
	double step_bis, t_bis, bissection_precision, bissection_step_min;
	double cross_condition, cross_sign;
	double surface_side_current, surface_side_previous;
	double surface_side_current_bis, surface_side_previous_bis;
	double m1_distance, collision_radius_big;
	double m2_distance, collision_radius_small;

	// initialize control variables
	bissection_precision = 1e-10;
	bissection_loop_max = 100;
	bissection_step_min = 1e-15;
	cross_condition = 0.0;
	cross_sign = 1.0;
	cross_coordinate = 1;
	cross_velocity = 3;
	map_coordinate = 0;
	map_velocity = 2;
	collision_radius_big = 1.66e-2;
	collision_radius_small = 4.52e-3;
	system_dimension = 4;

	// allocate memory and initializes exit data
	space_map = 1000;
	alloc_1d_double(energy, space_map);
	if (bissection_times != NULL)
	{
		alloc_1d_double(bissection_times, space_map);
	}
	alloc_2d_double(map, space_map, 2);

	// initialize map counter
	time_counter = 0;

	// FILE *tst = fopen("results/orbit/bissection_points.dat","w");

	// analyze every couple of adjoint orbit points
	for (int i = 1; i < step_number; i++)
	{
		// measures distance to primaries
		m1_distance = distance_to_primary(orbit[i], params, 0.0, "rotational", 1);
		m2_distance = distance_to_primary(orbit[i], params, 0.0, "rotational", 2);

		if (m1_distance < collision_radius_big ||
			m2_distance < collision_radius_small)
		{
			goto next;
		}

		// initialize cross status for both points
		surface_side_previous = orbit[i - 1][cross_coordinate] - cross_condition;
		surface_side_current = orbit[i][cross_coordinate] - cross_condition;

		// verify if points are on different sides of
		// the poincare surface
		if (surface_side_current * surface_side_previous < 0.0)
		{
			// initialize bissection parameters
			bissection_loop_control = 0;
			copy(y_bis, orbit[i - 1], system_dimension);
			t_bis = time[i - 1];
			step_bis = time[i] - time[i - 1];

		// return line if cross
		return_bissection_cross:;

			// update step size
			step_bis /= 2.0;
			if (fabs(step_bis) < bissection_step_min)
			{
				printf("Warning: minimun step reached at bissection.\n");
				break;
			}

		// return line if no cross
		return_bissection_no_cross:;

			// fprintf(tst, "%1.15e %1.15e\n", y_bis[0], y_bis[1]);

			// update bissection loop control
			bissection_loop_control++;
			if (bissection_loop_control == bissection_loop_max)
			{
				printf("Warning: bissection loop exceeded.\n");
				exit(2);
			}

			// store orbit before evolving
			copy(y_bis_temp, y_bis, system_dimension);

			// calculates time evolution
			evolve(y_bis, params, step_bis, NULL, NULL, NULL, NULL);

			// update surface sides
			surface_side_previous_bis = y_bis_temp[cross_coordinate] - cross_condition;		
			surface_side_current_bis = y_bis[cross_coordinate] - cross_condition;
			t_bis += step_bis;

			// check cross
			if (surface_side_current_bis * surface_side_previous_bis < 0.0)
			{
				// check if desired precision is reached
				if (fabs(surface_side_current_bis) < bissection_precision)
				{
					goto exit_loop_bissection;
				}
				else
				{
					copy(y_bis, y_bis_temp, system_dimension);
					t_bis -= step_bis;
					goto return_bissection_cross;
				}
			}
			else
			{
				goto return_bissection_no_cross;
			}

		// go to line if calculation was completed 
		exit_loop_bissection:;

			// check if cross is in the right direction
			if (cross_sign * y_bis[cross_velocity] > 0.0)
			{
				// allocate more memory space to exit data if necessary
				if (time_counter == space_map)
				{
					space_map += 1000;
					realloc_1d_double(energy, space_map);
					if (bissection_times != NULL)
					{
						realloc_1d_double(bissection_times, space_map);
					}
					realloc_2d_double(map, space_map, 2, time_counter);
				}

				// update exit data
				(*energy)[time_counter] = motion_integral(y_bis, params, t_bis, "rotational");
				if (bissection_times != NULL)
				{
					(*bissection_times)[time_counter] = t_bis;
				}
				(*map)[time_counter][0] = y_bis[map_coordinate];
				(*map)[time_counter][1] = y_bis[map_velocity];

				// update counter
				time_counter++;
			}

			// fprintf(tst, "\n");
		}

		next:;

	}

	// define exit data size
	*step_number_map = time_counter;
	*size_map = space_map;

	// fclose(tst);

	return 0;
}

int field_zero_velociy_curves_rotational(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
	f[0] = 2. * omega_y;
	f[1] = -2. * omega_x;
	return GSL_SUCCESS;
}

double root_function_L3(double x, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	double J = p->J;
	return x * x + 2. * (1. - mu) / fabs(x + mu) + 2. * mu / fabs(x + mu - 1.) - J;
}

double root_derivative_L3(double x, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	return 2. * x - 2. * (1. - mu) * (x + mu) / (fabs(x + mu) * fabs(x + mu) * fabs(x + mu)) - 2. * mu * (x + mu - 1.) / (fabs(x + mu - 1.) * fabs(x + mu - 1.) * fabs(x + mu - 1.));
}

void root_fdf_L3(double x, void *params, double *y, double *dy)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	double J = p->J;
	*y = x * x + 2. * (1. - mu) / fabs(x + mu) + 2. * mu / fabs(x + mu - 1.) - J;
	*dy = 2. * x - 2. * (1. - mu) * (x + mu) / (fabs(x + mu) * fabs(x + mu) * fabs(x + mu)) - 2. * mu * (x + mu - 1.) / (fabs(x + mu - 1.) * fabs(x + mu - 1.) * fabs(x + mu - 1.));
}

double root_function_triangular_points(double y, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	double J = p->J;
	double pos[4];
	pos[0] = 0.5 - mu;
	pos[1] = y;
	pos[2] = 0.0;
	pos[3] = 0.0;
	return 2. * pseudo_potential(pos, &mu) - J; 
}

double root_derivative_triangular_points(double y, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	double pos[4];
	pos[0] = 0.5 - mu;
	pos[1] = y;
	pos[2] = 0.0;
	pos[3] = 0.0;
	return 2. * pseudo_potential_y_derivative(pos, &mu); 
}

void root_fdf_triangular_points(double y, void *params, double *z, double *dz)
{
	struct root_params *p = (struct root_params *)params;
	double mu = p->mu;
	double J = p->J;
	double pos[4];
	pos[0] = 0.5 - mu;
	pos[1] = y;
	pos[2] = 0.0;
	pos[3] = 0.0;
	*z = 2. * pseudo_potential(pos, &mu) - J; 
	*dz = 2. * pseudo_potential_y_derivative(pos, &mu);
}

int root_zvc(double J, void *params, double *y, double x)
{
	int iter, max_iter, status;
	double mu, x0;
	mu = *(double *)params;
	struct root_params params_root = {mu, J};
	const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
	gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
	gsl_function_fdf FDF;
	if (J >= lagrangian_point_jacobi_constant(&mu, 3))
	{
		FDF.f = &root_function_L3;
		FDF.df = &root_derivative_L3;
		FDF.fdf = &root_fdf_L3;
	}
	else if (J >= lagrangian_point_jacobi_constant(&mu, 4))
	{
		FDF.f = &root_function_triangular_points;
		FDF.df = &root_derivative_triangular_points;
		FDF.fdf = &root_fdf_triangular_points;
	}
	FDF.params = &params_root;
	max_iter = 100;
	iter = 0;
	gsl_root_fdfsolver_set(s, &FDF, x);
	do
	{
		iter++;
		gsl_root_fdfsolver_iterate(s);
		x0 = x;
		x = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta(x, x0, 1e-15, 0);
		if (iter == max_iter)
		{
			printf("Warning: maximum iterate number reached at ZVC root\n");
			return 1;
		}
	} while (status == GSL_CONTINUE);
	if (J >= lagrangian_point_jacobi_constant(&mu, 3))
	{
		y[0] = x;
		y[1] = 0.;
	}
	else if (J >= lagrangian_point_jacobi_constant(&mu, 4))
	{
		y[0] = 0.5 - mu;
		y[1] = x;
	}
	gsl_root_fdfsolver_free(s);
	return 0;
}

int find_periodic_orbit_initial_parameters(void *params, double *y0, double half_period, double *constant, double *period)
{
	int count;
	double mu = *(double *)params, mu_1, mu_2, y[20];
	mu_1 = 1.0 - mu;
	mu_2 = mu;
	count = 0;
	while (count < 100)
	{
		copy(y, y0, 4);
		identity_matrix_array_form(y, 4, 4);
		evolve_extended(y, params, half_period);
		double r_1 =
		sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
		double r_2 =
		sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
		double omega_x =
		y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
		double g_2 = y[3];
		double g_3 = omega_x + 2. * y[3];
		double det_trans_matrix = y[15] * g_2 - y[11] * g_3;
		y0[3] -= (g_2 * y[2] - g_3 * y[1]) / det_trans_matrix;
		half_period -= (y[15] * y[1] - y[11] * y[2]) / det_trans_matrix;
		count++;
	}
	*constant = motion_integral(y0, params, 0.0, "extended");
	*period = 2.0 * half_period;
	// printf("J_upo = %1.5e period = %1.5e\n", *constant, *period);
	return 0;
}

int Lyapunov_orbit_initial_condition_linear_region(void *params, double *y0, double delta, double *J, double *period, int n)
{
	double mu = *(double *)params, gm, c2, om, kp, half_period;
	if(n == 1)
	{
		double L1x =
		collinear_lagrangian_point_x_coordinate(params, 1);
		gm = fabs(1. - mu - L1x);
		c2 = (mu + (1. - mu) * (gm / (1. - gm)) * (gm / (1. - gm)) * (gm / (1. - gm))) / (gm * gm * gm);
		om = sqrt(.5 * (2. - c2 + sqrt(9. * c2 * c2 - 8. * c2)));
		kp = (om * om + 2. * c2 + 1.) / (2. * om);
		y0[0] = L1x + delta;
		y0[1] = 0.;
		y0[2] = 0.;
		y0[3] = -kp * om * (y0[0] - (1. - mu) + gm);
	}
	else if(n == 2)
	{
		double L2x =
		collinear_lagrangian_point_x_coordinate(params, 2);
		gm = fabs(1. - mu - L2x);
		c2 = (mu + (1. - mu) * (gm / (1. + gm)) * (gm / (1. + gm)) * (gm / (1. + gm))) / (gm * gm * gm);
		om = sqrt(.5 * (2. - c2 + sqrt(9. * c2 * c2 - 8. * c2)));
		kp = (om * om + 2. * c2 + 1.) / (2. * om);
		y0[0] = L2x + delta;
		y0[1] = 0.;
		y0[2] = 0.;
		y0[3] = -kp * om * (y0[0] - (1. - mu) - gm);
	}
	else
	{
		printf("Warning:function not yet implemented for L3, L4, L5\n");
		exit(2);
	}	
	half_period = 1.35;
	find_periodic_orbit_initial_parameters(params, y0, half_period, J, period);
	return 0;
}

int find_Lyapunov_orbit_parameters(void *params, double Jupo, double *y0upo, double *periodupo, int n)
{
	// declare variables
	double constant_precision;
	double y0[4], y01[4], y02[4], yk[4], y03;
	double half_period_try, J, period;
	double delta_x, step_x, periodk;
	double delta_ydot, delta_ydotk;

	// define method parameters
	constant_precision = 1e-13;
	delta_x = 1e-5;
	step_x = 10. * delta_x;

	// calculate reference orbits in linear region
	Lyapunov_orbit_initial_condition_linear_region(params, y01, delta_x, &J, &period, n);
	Lyapunov_orbit_initial_condition_linear_region(params, y02, step_x, &J, &period, n);

	// check if motion constant value is adequate
	if (J < Jupo)
	{
		printf("Warning: There is no Lyapunov orbit for the given constant\n");
		exit(2);
	}
	else
	{
		printf("Searching for the Lyapunov orbit for the given constant\n");
	}

	// start method
	delta_ydot = fabs(y02[3] - y01[3]);
	copy(y0, y02, 4);
	do
	{
		do
		{
			y03 = y0[3];
			copy(yk, y0, 4);
			delta_ydotk = delta_ydot;
			periodk = period;
			y0[0] += step_x;
			y0[1] = 0.0;
			y0[2] = 0.0;
			y0[3] -= delta_ydot;
			half_period_try = period / 2.0;
			find_periodic_orbit_initial_parameters(params, y0, half_period_try, &J, &period);
			delta_ydot = fabs(y0[3] - y03);
		} while (J > Jupo);
		copy(y0, yk, 4);
		period = periodk;
		step_x /= 2.0;
		delta_ydot = delta_ydotk / 2.0;
		printf("Constant precision = %1.15e\n", fabs(J - Jupo));
	} while (fabs(J - Jupo) > constant_precision);

	// update initial condition and period of orbit
	copy(y0upo, y0, 4);
	*periodupo = period;

	return 0;
}

double determinant(double data[])
{
	// declare variables
	int s;
	double local_data[16], det;

	// make a local copy
	copy(local_data, data, 16);

	// transform no matrix representation
	gsl_matrix_view m =
	gsl_matrix_view_array(local_data, 4, 4);

	// calculate matrix decomposition
	gsl_permutation *p = gsl_permutation_alloc(4);
	gsl_linalg_LU_decomp(&m.matrix, p, &s);
	
	// calculate determinant
	det = gsl_linalg_LU_det(&m.matrix, s);
	
	// free memory
	gsl_permutation_free(p);
	
	return det;
}

int calculate_eigen_values_and_vectors(double data[], double eig_val_real[], double eig_val_imag[], double eig_vec_real[][4], double eig_vec_imag[][4])
{
	// declare variables
	double local_data[16];

	// make a local copy
	copy(local_data, data, 16);

	// prints matrix determinant
	// printf("Determinant = %1.5e\n", determinant(local_data));
	
	// creates matrix representation
	gsl_matrix_view m = gsl_matrix_view_array(local_data, 4, 4);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(4);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(4, 4);
	
	// define workspace
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(4);
	
	// calculate eigenvalues and eigenvectors
	gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);

	// sort eigenvectors by descending order of eigenvalues
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	// print results
	// for (int i = 0; i < 4; i++)
	// {
	// 	gsl_complex eval_i = gsl_vector_complex_get(eval, i);
	// 	gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
	// 	printf("eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
	// 	printf("eigenvector = \n");
	// 	for (int j = 0; j < 4; ++j)
	// 	{
	// 		gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	// 		printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
	// 	}
	// }

	// get eigenvalues and eigenvectors
	// and pass results
	for (int i = 0; i < 4; i++)
	{
		gsl_complex eval_i = gsl_vector_complex_get(eval, i);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

		eig_val_real[i] = GSL_REAL(eval_i);
		eig_val_imag[i] = GSL_IMAG(eval_i);
		for (int j = 0; j < 4; ++j)
		{
			gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
			eig_vec_real[i][j] = GSL_REAL(z);
			eig_vec_imag[i][j] = GSL_IMAG(z);
		}
	}

	// free memory
	gsl_eigen_nonsymmv_free(w);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);

	return 0;
}

int find_outer_lunar_orbit(void *params, double J_orbit, double *y0_orbit, double *period_orbit, char* stability)
{
	// declare variables
	double constant_precision;
	double y0[4], y01[4], y02[4], yk[4];
	double y03, J0, dir;
	double half_period_try, J, period;
	double step_x, periodk;
	double delta_ydot, delta_ydotk;

	// define method parameters
	constant_precision = 1e-12; // 1e-13

	// check if motion constant value is adequate
	if (J_orbit > 3.18266)
	{
		printf("Warning: There is no outer lunar stable orbit for the given constant\n");
		exit(2);
	}

	if (strcmp(stability, "stable") == 0)
	{
		// define step for method
		step_x = 1e-3;

		// define direction
		dir = 1.0;

		// calculate reference orbits
		y01[0] = 1.11089;
		y01[1] = 0.;
		y01[2] = 0.;
		integral_to_velocity(y01, params, 0.0, 3.1730020406, "rotational", 3);
		find_periodic_orbit_initial_parameters(params, y01, 0.956, &J0, &period);
		y02[0] = 1.11079;
		y02[1] = 0.;
		y02[2] = 0.;
		integral_to_velocity(y02, params, 0.0, 3.1730225223, "rotational", 3);
		find_periodic_orbit_initial_parameters(params, y02, 0.956, &J0, &period);
	}
	else if (strcmp(stability, "unstable") == 0)
	{
		// define direction
		dir = -1.0;

		// define step for method
		if (J_orbit > 3.1789)
		{
			step_x = 1e-4;
		}
		else if (J_orbit > 3.1765)
		{
			step_x = 1e-6;
		}
		else
		{
			step_x = 1e-6;
			printf("function not working for this J\n");
			exit(2);
		}
		
		// calculate reference orbits
		y02[0] = 1.0607;
		y02[1] = 0.;
		y02[2] = 0.;
		integral_to_velocity(y02, params, 0.0, 3.17951, "rotational", 3);
		find_periodic_orbit_initial_parameters(params, y02, 0.703945, &J0, &period);

		y01[0] = 1.0606;
		y01[1] = 0.;
		y01[2] = 0.;
		integral_to_velocity(y01, params, 0.0, 3.17925, "rotational", 3);
		find_periodic_orbit_initial_parameters(params, y01, 0.703945, &J0, &period);

	}
	else
	{
		printf("Warning: wrong stability for outer lunar orbit calculation\n");
		exit(2);
	}
	
	// start method
	delta_ydot = fabs(y02[3] - y01[3]);
	copy(y0, y01, 4);
	do
	{
		do
		{
			y03 = y0[3];
			copy(yk, y0, 4);
			delta_ydotk = delta_ydot;
			periodk = period;
			y0[0] += dir * sign(J0 - J_orbit) * step_x;
			y0[1] = 0.0;
			y0[2] = 0.0;
			y0[3] -= dir * sign(J0 - J_orbit) * delta_ydot;
			half_period_try = period / 2.0;
			find_periodic_orbit_initial_parameters(params, y0, half_period_try, &J, &period);
			delta_ydot = fabs(y0[3] - y03);
		} while (sign(J0 - J_orbit) * (J_orbit - J) < 0.0);
		copy(y0, yk, 4);
		period = periodk;
		step_x /= 2.0;
		delta_ydot = delta_ydotk / 2.0;
		printf("Orbit precision = %1.15e\n", fabs(J - J_orbit));
	} while (fabs(J - J_orbit) > constant_precision);

	// update initial condition and period of orbit
	copy(y0_orbit, y0, 4);
	*period_orbit = period;

	// printf("%1.3f %1.10e\n", J_orbit, period);

	return 0;
}

int profile_outer_lunar_orbit(void *params)
{
	FILE *spos = fopen("results/lunar_periodic_orbits/outer_orbit_stable_pos.dat", "w");
	FILE *upos = fopen("results/lunar_periodic_orbits/outer_orbit_unstable_pos.dat", "w");
	FILE *seig = fopen("results/lunar_periodic_orbits/outer_orbit_stable_eig.dat", "w");
	FILE *ueig0 = fopen("results/lunar_periodic_orbits/outer_orbit_unstable_eig0.dat", "w");
	FILE *ueig1 = fopen("results/lunar_periodic_orbits/outer_orbit_unstable_eig1.dat", "w");

	FILE *spos_J, *upos_J;
	char filename[100];

	int count;
	double J, J_max, J_min, delta_J;
	double y0[4], y[20], period;
	double monodromy_matrix[16];
	double eig_val_real[4], eig_val_imag[4];
	double eig_vec_real[4][4], eig_vec_imag[4][4];
	double eig_val[2], eig_vec[2][4];

	// J_max = 3.1826604103108;
	// J_min = 3.182660;
	// delta_J = 1e-7;

	J_max = 3.176;
	J_min = 3.1755;
	delta_J = 1e-3;

	J = J_max;
	while (J > J_min)
	{
		printf("J = %1.5f\n", J);

		printf("Calculating stable orbit\n");

		find_outer_lunar_orbit(params, J, y0, &period, "stable");

		fprintf(spos, "%1.15e %1.15e\n", J, y0[0]);

		sprintf(filename, "results/lunar_periodic_orbits/outer_orbit_stable_pos_%1.3f.dat", J);	
		spos_J = fopen(filename, "w");
		fprintf(spos_J, "%1.15e %1.15e\n", y0[0], 0.0);
		fclose(spos_J);

		copy(y, y0, 4);
		identity_matrix_array_form(y, 4, 4);
		evolve_extended(y, params, period);

		for (int i = 0; i < 16; i++)
		{
			monodromy_matrix[i] = y[i + 4];
		}

		calculate_eigen_values_and_vectors(monodromy_matrix, eig_val_real, eig_val_imag, eig_vec_real, eig_vec_imag);

		count = 0;
		for (int i = 0; i < 4; i++)
		{
			if (fabs(eig_val_real[i] - 1.0) > 1e-3)
			{
				eig_val[count] = eig_val_real[i];
				copy(eig_vec[count], eig_vec_real[i], 4);

				count++;
			}
		}

		fprintf(seig, "%1.15e %1.15e\n", J, eig_val[0]);
		fprintf(seig, "%1.15e %1.15e\n", J, eig_val[1]);

		printf("Calculating unstable orbit\n");

		find_outer_lunar_orbit(params, J, y0, &period, "unstable");

		fprintf(upos, "%1.15e %1.15e\n", J, y0[0]);

		sprintf(filename, "results/lunar_periodic_orbits/outer_orbit_unstable_pos_%1.3f.dat", J);	
		upos_J = fopen(filename, "w");
		fprintf(upos_J, "%1.15e %1.15e\n", y0[0], 0.0);
		fclose(upos_J);

		copy(y, y0, 4);
		identity_matrix_array_form(y, 4, 4);
		evolve_extended(y, params, period);

		for (int i = 0; i < 16; i++)
		{
			monodromy_matrix[i] = y[i + 4];
		}

		calculate_eigen_values_and_vectors(monodromy_matrix, eig_val_real, eig_val_imag, eig_vec_real, eig_vec_imag);

		count = 0;
		for (int i = 0; i < 4; i++)
		{
			// separate those with non-unity real part
			if (fabs(eig_val_real[i] - 1.0) > 1e-2)
			{

				// pass results
				eig_val[count] = eig_val_real[i];
				copy(eig_vec[count], eig_vec_real[i], 4);

				// change result array index
				count++;
			}
		}

		fprintf(ueig0, "%1.15e %1.15e\n", J, eig_val[0]);
		fprintf(ueig1, "%1.15e %1.15e\n", J, eig_val[1]);

		J -= delta_J;
	}

	fclose(spos);
	fclose(upos);
	fclose(seig);
	fclose(ueig0);
	fclose(ueig1);

	return 0;
}

int profile_outer_lunar_orbit_v2(void *params)
{
	FILE *upos = fopen("results/lunar_periodic_orbits/outer_orbit_unstable_pos.dat", "w");

	double y0[4];
	double J_orbit, period;
	double x_max, x_min, delta_x;

	x_min = 1.06033;
	x_max = 1.06250;
	delta_x = 1e-5;

	y0[0] = x_max;
	y0[1] = 0.0;
	y0[2] = 0.0;
	integral_to_velocity(y0, params, 0.0, 3.186, "rotational", 3);

	find_periodic_orbit_initial_parameters(params, y0, 0.703945, &J_orbit, &period);

	for (double x = x_max; x > x_min; x -= delta_x) 
	{
		printf("Calculating orbit at x = %1.8f\n", x);

		if (x < 1.0604)
		{
			delta_x = 5e-6;
		}

		if (x < 1.06036)
		{
			delta_x = 1e-6;
		}

		// if (x < 1.0603509)
		// {
		// 	delta_x = 1e-8;
		// }

		y0[0] = x;
		y0[1] = 0.0;
		y0[2] = 0.0;
		integral_to_velocity(y0, params, 0.0, J_orbit, "rotational", 3);

		find_periodic_orbit_initial_parameters(params, y0, period / 2.0, &J_orbit, &period);

		fprintf(upos, "%1.15e %1.15e\n", J_orbit, y0[0]);
	}

	fclose(upos);
	
	return 0;
}

int profile_inner_lunar_orbit(void *params)
{
	FILE *pos = fopen("results/lunar_periodic_orbits/inner_orbit_pos.dat", "w");
	FILE *eig0 = fopen("results/lunar_periodic_orbits/inner_orbit_eig0.dat", "w");
	FILE *eig1 = fopen("results/lunar_periodic_orbits/inner_orbit_eig1.dat", "w");

	int count;
	double y0[4], y[20];
	double J_orbit, period;
	double monodromy_matrix[16];
	double eig_val[2], eig_vec[2][4];
	double eig_val_real[4], eig_val_imag[4];
	double eig_vec_real[4][4], eig_vec_imag[4][4];
	double x_max, x_min, delta_x;

	x_min = 1.000;
	x_max = 1.016;
	delta_x = 5e-5;

	for (double x = x_max; x > x_min; x -= delta_x) 
	{
		printf("Calculating orbit at x = %1.5f\n", x);

		y0[0] = x;
		y0[1] = 0.0;
		y0[2] = 0.0;
		integral_to_velocity(y0, params, 0.0, 3.186, "rotational", 3);

		find_periodic_orbit_initial_parameters(params, y0,  0.864, &J_orbit, &period);

		fprintf(pos, "%1.15e %1.15e\n", J_orbit, y0[0]);

		copy(y, y0, 4);
		identity_matrix_array_form(y, 4, 4);
		evolve_extended(y, params, period);

		for (int i = 0; i < 16; i++)
		{
			monodromy_matrix[i] = y[i + 4];
		}

		calculate_eigen_values_and_vectors(monodromy_matrix, eig_val_real, eig_val_imag, eig_vec_real, eig_vec_imag);

		count = 0;
		for (int i = 0; i < 4; i++)
		{
			if (fabs(eig_val_real[i] - 1.0) > 1e-1)
			{
				eig_val[count] = eig_val_real[i];
				copy(eig_vec[count], eig_vec_real[i], 4);

				count++;
			}
		}

		fprintf(eig0, "%1.15e %1.15e\n", J_orbit, eig_val[0]);
		fprintf(eig1, "%1.15e %1.15e\n", J_orbit, eig_val[1]);
	}

	fclose(pos);
	fclose(eig0);
	fclose(eig1);
	
	return 0;
}

int profile_inner_lunar_orbit_v2(void *params)
{
	FILE *pos = fopen("results/lunar_periodic_orbits/inner_orbit_pos.dat", "w");

	double y0[4];
	double J_orbit, period;
	double x_max, x_min, delta_x;

	x_min = 0.99305;
	x_max = 0.9932;
	delta_x = 1e-5;

	for (double x = x_max; x > x_min; x -= delta_x) 
	{
		printf("Calculating orbit at x = %1.5f\n", x);

		y0[0] = x;
		y0[1] = 0.0;
		y0[2] = 0.0;
		integral_to_velocity(y0, params, 0.0, 3.186, "rotational", 3);

		find_periodic_orbit_initial_parameters(params, y0,  0.864, &J_orbit, &period);

		fprintf(pos, "%1.15e %1.15e\n", J_orbit, y0[0]);

	}

	fclose(pos);

	return 0;
}

int trace_manifolds_outer_saddle(double J_orbit)
{
	// declare exit files
	FILE *mupm, *mumm, *mspm, *msmm;
	FILE *mupmj, *mummj;
	char filename[100];

	// open exit files

	// manifolds on poincare map
	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_unstable_plus.dat");
	mupm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_unstable_minus.dat");
	mumm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_stable_plus.dat");
	mspm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_stable_minus.dat");	
	msmm = fopen(filename, "w");

	// jacobi constant of 
	// manifolds on poincare map
	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_unstable_plus_jacobi.dat");
	mupmj = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_outer_saddle_unstable_minus_jacobi.dat");
	mummj = fopen(filename, "w");

	// declare variables
	double mu = 1.215e-2;
	int man_ic_number, count;
	int size_orbit_mum, size_orbit_mup;
	int size_map_mumm, size_map_mupm;
	int step_number_mum, step_number_mup;
	int step_number_mumm, step_number_mupm;
	double y0[4], y[20], **man_ic;
	double man_ic_dist_from_upo, t_man_minus, t_man_plus;
	double period, man_ic_step, norm;
	double monodromy_matrix[16], eig_val[2], eig_vec[2][4];
	double eig_val_real[4], eig_val_imag[4];
	double eig_vec_real[4][4], eig_vec_imag[4][4];
	double transition_matrix[16], eig_vec_uns_evolve[4];
	double man_ic_uns_plus[4], man_ic_uns_minus[4];
	double **orbit_mum, *time_mum, **map_mum, *map_energy_mum;
	double **orbit_mup, *time_mup, **map_mup, *map_energy_mup;

	// define function parameters
	man_ic_dist_from_upo = 1e-6;
	t_man_minus = 75.;
	t_man_plus = 32.;
	man_ic_number = 10000;

	// allocate memory
	alloc_2d_double(&man_ic, man_ic_number, 20);

	// calculate period-7 upo
	find_outer_lunar_orbit(&mu, J_orbit, y0, &period, "unstable");
	copy(y, y0, 4);
	identity_matrix_array_form(y, 4, 4);

	printf("Here\n");

	// determine distance between 
	// points of Lyapunov orbit
	man_ic_step = period / (double)man_ic_number;

	// calculate points on Lyapunov orbit
	// along the transition matrices
	for (int i = 0; i < man_ic_number; i++)
	{
		evolve_extended(y, &mu, man_ic_step);
		copy(man_ic[i], y, 20);
	}

	// define monodromy matrix
	for (int i = 0; i < 16; i++)
	{
		monodromy_matrix[i] = y[i + 4];
	}

	// calculate eigenvalues and eigenvectors
	// of the monodromy matrix
	calculate_eigen_values_and_vectors(monodromy_matrix, eig_val_real, eig_val_imag, eig_vec_real, eig_vec_imag);

	count = 0;
	for (int i = 0; i < 4; i++)
	{
		// separate those with non-unity real part
		if (fabs(eig_val_real[i] - 1.0) > 1e-2)
		{
			// separate those with no imaginary part
			if(fabs(eig_val_imag[i]) < 1e-2)
			{
				// pass results
				eig_val[count] = eig_val_real[i];
				copy(eig_vec[count], eig_vec_real[i], 4);

				// change result array index
				count++;
			}
		}
	}

	// check eigenvalues
	if (eig_val[0] < 1.0 || eig_val[1] > 1.0)
	{
		printf("Warning: something wrong with eigenvalues calculation in Lyapunov orbit manifolds tracing\n");
		printf("%d %1.3e %1.3e\n", count, eig_val[0], eig_val[1]);
		print_array(eig_val_real,4);
		print_array(eig_val_imag,4);
		exit(2);
	}

	// loop on Lyapunov orbit points
	for (int i = 0; i < man_ic_number; i++)
	{
		// print progress
		printf("Calculating orbit %d of %d\n", i + 1, man_ic_number);

		// define transition matrix
		for (int j = 0; j < 16; j++)
		{
			transition_matrix[j] = man_ic[i][j + 4];
		}

		// calculate point eigenvectors using
		// the monodromy matrix eigenvectors
		// and the point transition matrix
		gsl_matrix_const_view TM =
		gsl_matrix_const_view_array(transition_matrix, 4, 4);
		gsl_vector_const_view EVU =
		gsl_vector_const_view_array(eig_vec[0], 4);
		gsl_vector_view EVUE =
		gsl_vector_view_array(eig_vec_uns_evolve, 4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &TM.matrix, &EVU.vector, 0.0, &EVUE.vector);

		// normalize point eigenvectors
		norm = array_norm(eig_vec_uns_evolve, 4);
		for (int j = 0; j < 4; j++)
		{
			eig_vec_uns_evolve[j] /= norm;
		}

		// calculate manifold initial condition 
		for (int j = 0; j < 4; j++)
		{
			man_ic_uns_plus[j] = man_ic[i][j] + man_ic_dist_from_upo * eig_vec_uns_evolve[j];
			man_ic_uns_minus[j] = man_ic[i][j] - man_ic_dist_from_upo * eig_vec_uns_evolve[j];
		}

		// calculate manifold orbit minus branch
		evolve(man_ic_uns_minus, &mu, t_man_minus, &orbit_mum, &time_mum, &step_number_mum, &size_orbit_mum);

		// calculate manifold orbit on poincare map
		poincare_map(&mu, &map_mum, &map_energy_mum, NULL, &step_number_mumm, &size_map_mumm, orbit_mum, time_mum, step_number_mum);
		
		// write manifold orbit minus branch on poincare map
		for (int j = 0; j < step_number_mumm; j++)
		{
			// write unstable manifold on poincare map
			fprintf(mumm, "%1.15e %1.15e\n", map_mum[j][0], map_mum[j][1]);

			// write stable manifold on poincare map
			fprintf(msmm, "%1.15e %1.15e\n", map_mum[j][0], -map_mum[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mummj, "%d %1.15e\n", j, fabs(map_energy_mum[j] - J_orbit));
		}
		fprintf(mumm, "\n");
		fprintf(msmm, "\n");
		fprintf(mummj, "\n");

		// free memory minus branch
		dealloc_1d_double(&time_mum);
		dealloc_2d_double(&orbit_mum, size_orbit_mum);
		dealloc_1d_double(&map_energy_mum);
		dealloc_2d_double(&map_mum, size_map_mumm);

		// calculate manifold orbit plus branch
		evolve(man_ic_uns_plus, &mu, t_man_plus, &orbit_mup, &time_mup, &step_number_mup, &size_orbit_mup);

		// calculate manifold orbit on poincare map
		poincare_map(&mu, &map_mup, &map_energy_mup, NULL, &step_number_mupm, &size_map_mupm, orbit_mup, time_mup, step_number_mup);

		// write manifold orbit plus branch on poincare map
		for (int j = 0; j < step_number_mupm; j++)
		{
			// write unstable manifold on poincare map
			fprintf(mupm, "%1.15e %1.15e\n", map_mup[j][0], map_mup[j][1]);

			// write stable manifold on poincare map
			fprintf(mspm, "%1.15e %1.15e\n", map_mup[j][0], -map_mup[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mupmj, "%d %1.15e\n", j, fabs(map_energy_mup[j] - J_orbit));
		}
		fprintf(mupm, "\n");
		fprintf(mspm, "\n");
		fprintf(mupmj, "\n");

		// free memory plus branch
		dealloc_1d_double(&time_mup);
		dealloc_2d_double(&orbit_mup, size_orbit_mup);
		dealloc_1d_double(&map_energy_mup);
		dealloc_2d_double(&map_mup, size_map_mupm);
	}

	// free memory
	dealloc_2d_double(&man_ic, man_ic_number);

	// close exit files
	fclose(mupm); fclose(mumm); fclose(mspm); fclose(msmm);
	fclose(mupmj); fclose(mummj);

	return 0;
}

int find_satellite_UPO(double y0[4], double *J, double *period, char *island)
{
	double mu = 1.215e-2;
	double J_island, half_period_try;

	if (strcmp(island, "inner") == 0)
	{
		J_island = 3.187;
		half_period_try = 5.345;
		y0[0] = 1.0240819696;
	}
	else if (strcmp(island, "outer") == 0)
	{
		J_island = 3.176;
		half_period_try = 12.6155 / 2.0;
		y0[0] = 1.069955;
	}
	else
	{
		printf("Warning: wrong island for period 7 UPO determination\n");
		exit(2);
	}

	y0[1] = 0.0;
	y0[2] = 0.0;
	integral_to_velocity(y0, &mu, 0.0, J_island, "rotational", 3);

	find_periodic_orbit_initial_parameters(&mu, y0, half_period_try, J, period);

	return 0;
}

int trace_manifolds_satellite_UPO(double *J_orbit, char *island)
{
	// declare exit files
	FILE *mupm, *mumm, *mspm, *msmm;
	FILE *mupmj, *mummj;
	char filename[100];

	// open exit files

	// manifolds on poincare map
	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_unstable_plus_%s.dat", island);
	mupm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_unstable_minus_%s.dat", island);
	mumm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_stable_plus_%s.dat", island);
	mspm = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_stable_minus_%s.dat", island);	
	msmm = fopen(filename, "w");

	// jacobi constant of 
	// manifolds on poincare map
	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_unstable_plus_jacobi_%s.dat", island);
	mupmj = fopen(filename, "w");

	sprintf(filename, "results/lunar_periodic_orbits/manifold_satellite_UPO_unstable_minus_jacobi_%s.dat", island);
	mummj = fopen(filename, "w");

	// declare variables
	double mu = 1.215e-2;
	int man_ic_number, count;
	int size_orbit_mum, size_orbit_mup;
	int size_map_mumm, size_map_mupm;
	int step_number_mum, step_number_mup;
	int step_number_mumm, step_number_mupm;
	double y0[4], y[20], **man_ic;
	double man_ic_dist_from_upo, t_man_minus, t_man_plus;
	double period, man_ic_step, norm;
	double monodromy_matrix[16], eig_val[2], eig_vec[2][4];
	double eig_val_real[4], eig_val_imag[4];
	double eig_vec_real[4][4], eig_vec_imag[4][4];
	double transition_matrix[16], eig_vec_uns_evolve[4];
	double man_ic_uns_plus[4], man_ic_uns_minus[4];
	double **orbit_mum, *time_mum, **map_mum, *map_energy_mum;
	double **orbit_mup, *time_mup, **map_mup, *map_energy_mup;

	// define function parameters
	man_ic_dist_from_upo = 1e-6;
	t_man_minus = 200.;
	t_man_plus = 200.;
	man_ic_number = 10000;

	// allocate memory
	alloc_2d_double(&man_ic, man_ic_number, 20);

	// calculate period-7 upo
	find_satellite_UPO(y0, J_orbit, &period, island);
	copy(y, y0, 4);
	identity_matrix_array_form(y, 4, 4);

	// determine distance between 
	// points of Lyapunov orbit
	man_ic_step = period / (double)man_ic_number;

	// calculate points on Lyapunov orbit
	// along the transition matrices
	for (int i = 0; i < man_ic_number; i++)
	{
		evolve_extended(y, &mu, man_ic_step);
		copy(man_ic[i], y, 20);
	}

	// define monodromy matrix
	for (int i = 0; i < 16; i++)
	{
		monodromy_matrix[i] = y[i + 4];
	}

	// calculate eigenvalues and eigenvectors
	// of the monodromy matrix
	calculate_eigen_values_and_vectors(monodromy_matrix, eig_val_real, eig_val_imag, eig_vec_real, eig_vec_imag);

	count = 0;
	for (int i = 0; i < 4; i++)
	{
		// separate those with non-unity real part
		if (fabs(eig_val_real[i] - 1.0) > 1e-2)
		{
			// separate those with no imaginary part
			if(fabs(eig_val_imag[i]) < 1e-2)
			{
				// pass results
				eig_val[count] = eig_val_real[i];
				copy(eig_vec[count], eig_vec_real[i], 4);

				// change result array index
				count++;
			}
		}
	}

	// check eigenvalues
	if (eig_val[0] < 1.0 || eig_val[1] > 1.0)
	{
		printf("Warning: something wrong with eigenvalues calculation in Lyapunov orbit manifolds tracing\n");
		printf("%d %1.3e %1.3e\n", count, eig_val[0], eig_val[1]);
		print_array(eig_val_real,4);
		print_array(eig_val_imag,4);
		exit(2);
	}

	// loop on Lyapunov orbit points
	for (int i = 0; i < man_ic_number; i++)
	{
		// print progress
		printf("Calculating orbit %d of %d\n", i + 1, man_ic_number);

		// define transition matrix
		for (int j = 0; j < 16; j++)
		{
			transition_matrix[j] = man_ic[i][j + 4];
		}

		// calculate point eigenvectors using
		// the monodromy matrix eigenvectors
		// and the point transition matrix
		gsl_matrix_const_view TM =
		gsl_matrix_const_view_array(transition_matrix, 4, 4);
		gsl_vector_const_view EVU =
		gsl_vector_const_view_array(eig_vec[0], 4);
		gsl_vector_view EVUE =
		gsl_vector_view_array(eig_vec_uns_evolve, 4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &TM.matrix, &EVU.vector, 0.0, &EVUE.vector);

		// normalize point eigenvectors
		norm = array_norm(eig_vec_uns_evolve, 4);
		for (int j = 0; j < 4; j++)
		{
			eig_vec_uns_evolve[j] /= norm;
		}

		// calculate manifold initial condition 
		for (int j = 0; j < 4; j++)
		{
			man_ic_uns_plus[j] = man_ic[i][j] + man_ic_dist_from_upo * eig_vec_uns_evolve[j];
			man_ic_uns_minus[j] = man_ic[i][j] - man_ic_dist_from_upo * eig_vec_uns_evolve[j];
		}

		// calculate manifold orbit minus branch
		evolve(man_ic_uns_minus, &mu, t_man_minus, &orbit_mum, &time_mum, &step_number_mum, &size_orbit_mum);

		// calculate manifold orbit on poincare map
		poincare_map(&mu, &map_mum, &map_energy_mum, NULL, &step_number_mumm, &size_map_mumm, orbit_mum, time_mum, step_number_mum);
		
		// write manifold orbit minus branch on poincare map
		for (int j = 0; j < step_number_mumm; j++)
		{
			// write unstable manifold on poincare map
			fprintf(mumm, "%1.15e %1.15e\n", map_mum[j][0], map_mum[j][1]);

			// write stable manifold on poincare map
			fprintf(msmm, "%1.15e %1.15e\n", map_mum[j][0], -map_mum[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mummj, "%d %1.15e\n", j, fabs(map_energy_mum[j] - *J_orbit));
		}
		fprintf(mumm, "\n");
		fprintf(msmm, "\n");
		fprintf(mummj, "\n");

		// free memory minus branch
		dealloc_1d_double(&time_mum);
		dealloc_2d_double(&orbit_mum, size_orbit_mum);
		dealloc_1d_double(&map_energy_mum);
		dealloc_2d_double(&map_mum, size_map_mumm);

		// calculate manifold orbit plus branch
		evolve(man_ic_uns_plus, &mu, t_man_plus, &orbit_mup, &time_mup, &step_number_mup, &size_orbit_mup);

		// calculate manifold orbit on poincare map
		poincare_map(&mu, &map_mup, &map_energy_mup, NULL, &step_number_mupm, &size_map_mupm, orbit_mup, time_mup, step_number_mup);

		// write manifold orbit plus branch on poincare map
		for (int j = 0; j < step_number_mupm; j++)
		{
			// write unstable manifold on poincare map
			fprintf(mupm, "%1.15e %1.15e\n", map_mup[j][0], map_mup[j][1]);

			// write stable manifold on poincare map
			fprintf(mspm, "%1.15e %1.15e\n", map_mup[j][0], -map_mup[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mupmj, "%d %1.15e\n", j, fabs(map_energy_mup[j] - *J_orbit));
		}
		fprintf(mupm, "\n");
		fprintf(mspm, "\n");
		fprintf(mupmj, "\n");

		// free memory plus branch
		dealloc_1d_double(&time_mup);
		dealloc_2d_double(&orbit_mup, size_orbit_mup);
		dealloc_1d_double(&map_energy_mup);
		dealloc_2d_double(&map_mup, size_map_mupm);
	}

	// free memory
	dealloc_2d_double(&man_ic, man_ic_number);

	// close exit files
	fclose(mupm); fclose(mumm); fclose(mspm); fclose(msmm);
	fclose(mupmj); fclose(mummj);

	return 0;
}