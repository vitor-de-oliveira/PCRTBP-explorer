#include "dynamical_system.h"

int field_extended(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double K1 = (1. - mu) / (r_1 * r_1 * r_1);
	double K2 = mu / (r_2 * r_2 * r_2);
	double A = (y[0] + mu) * (y[0] + mu) / (r_1 * r_1);
	double B = (y[0] + mu - 1.) * (y[0] + mu - 1.) / (r_2 * r_2);
	double C = y[1] * y[1] / (r_1 * r_1);
	double D = y[1] * y[1] / (r_2 * r_2);
	double E1 = (1. - mu) * (y[0] + mu) / (r_1 * r_1 * r_1 * r_1 * r_1);
	double E2 = mu * (y[0] + mu - 1.) / (r_2 * r_2 * r_2 * r_2 * r_2);
	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
	double omega_xx = 1. - K1 * (1. - 3. * A) - K2 * (1. - 3. * B);
	double omega_yy = 1. - K1 * (1. - 3. * C) - K2 * (1. - 3. * D);
	double omega_xy = 3. * y[1] * (E1 + E2);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = omega_x + 2.0 * y[3];
	f[3] = omega_y - 2.0 * y[2];
	double Df[] = {0., 0., 1., 0.,
				   0., 0., 0., 1.,
				   omega_xx, omega_xy, 0., 2.,
				   omega_xy, omega_yy, -2., 0.};
	double Phi[] = {y[4], y[5], y[6], y[7],
					y[8], y[9], y[10], y[11],
					y[12], y[13], y[14], y[15],
					y[16], y[17], y[18], y[19]};
	double F[] = {0., 0., 0., 0.,
				  0., 0., 0., 0.,
				  0., 0., 0., 0.,
				  0., 0., 0., 0.};
	gsl_matrix_view MDf = gsl_matrix_view_array(Df, 4, 4);
	gsl_matrix_view MPhi = gsl_matrix_view_array(Phi, 4, 4);
	gsl_matrix_view MF = gsl_matrix_view_array(F, 4, 4);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., &MDf.matrix, &MPhi.matrix, 0., &MF.matrix);
	for (int i = 0; i < 16; i++)
		f[i + 4] = F[i];
	return GSL_SUCCESS;
}

int field_hamiltonian(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
	f[0] = y[2] + y[1];
	f[1] = y[3] - y[0];
	f[2] = omega_x + y[3] - y[0];
	f[3] = omega_y - y[2] - y[1];
	return GSL_SUCCESS;
}

int field_inertial(double t, const double y[], double f[], void *params)
{
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2 * cos(t)) * (y[0] + mu_2 * cos(t)) + (y[1] + mu_2 * sin(t)) * (y[1] + mu_2 * sin(t)));
	double r_2 = sqrt((y[0] - mu_1 * cos(t)) * (y[0] - mu_1 * cos(t)) + (y[1] - mu_1 * sin(t)) * (y[1] - mu_1 * sin(t)));
	double omega_x = - mu_1 * (y[0] + mu_2 * cos(t)) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1 * cos(t)) / (r_2 * r_2 * r_2);
	double omega_y = - mu_1 * (y[1] + mu_2 * sin(t)) / (r_1 * r_1 * r_1) - mu_2 * (y[1] - mu_1 * sin(t)) / (r_2 * r_2 * r_2);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = omega_x;
	f[3] = omega_y;
	return GSL_SUCCESS;
}

int field_regularized_bigger_mass(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double *par = (double *)params;
	double mu = par[0];
	double J = par[1];
	double C = J + mu * (1.0 - mu);
	double u_squared = y[0] * y[0];
	double v_squared = y[1] * y[1];
	double squared_sum = u_squared + v_squared;
	double squared_dif = u_squared - v_squared;
	double D1_u = 6.0 * y[0] * squared_sum * squared_sum;
	double D1_v = 6.0 * y[1] * squared_sum * squared_sum;
	double D2_u = -4.0 * y[0] * y[0] * y[0];
	double D2_v = 4.0 * y[1] * y[1] * y[1];
	double D3_u = 2.0 * y[0];
	double D3_v = 2.0 * y[1];
	double D4_u_P1 = 2.0 * y[0] / sqrt(1.0 + squared_sum * squared_sum - 2.0 * squared_dif);
	double D4_u_P2 = 2.0 * y[0] * squared_sum * (squared_sum - 1.0) / pow(1.0 + squared_sum * squared_sum - 2.0 * squared_dif, 1.5);
	double D4_v_P1 = 2.0 * y[1] / sqrt(1.0 + squared_sum * squared_sum - 2.0 * squared_dif);
	double D4_v_P2 = 2.0 * y[1] * squared_sum * (squared_sum + 1.0) / pow(1.0 + squared_sum * squared_sum - 2.0 * squared_dif, 1.5);
	double D4_u = D4_u_P1 - D4_u_P2;
	double D4_v = D4_v_P1 - D4_v_P2;
	double potential_u = 2.0 * (D1_u + 2.0 * mu * D2_u + (mu - C) * D3_u + 2.0 * mu * D4_u);
	double potential_v = 2.0 * (D1_v + 2.0 * mu * D2_v + (mu - C) * D3_v + 2.0 * mu * D4_v);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = potential_u + 8.0 * squared_sum * y[3];
	f[3] = potential_v - 8.0 * squared_sum * y[2];
	return GSL_SUCCESS;
}

int field_regularized_smaller_mass(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double *par = (double *)params;
	double mu = par[0];
	double J = par[1];
	double C = J + mu * (1.0 - mu);
	double u_squared = y[0] * y[0];
	double v_squared = y[1] * y[1];
	double squared_sum = u_squared + v_squared;
	double squared_dif = u_squared - v_squared;
	double D1_u = 6.0 * y[0] * squared_sum * squared_sum;
	double D1_v = 6.0 * y[1] * squared_sum * squared_sum;
	double D2_u = 4.0 * y[0] * y[0] * y[0];
	double D2_v = -4.0 * y[1] * y[1] * y[1];
	double D3_u = 2.0 * y[0];
	double D3_v = 2.0 * y[1];
	double D4_u_P1 = 2.0 * y[0] / sqrt(1.0 + squared_sum * squared_sum + 2.0 * squared_dif);
	double D4_u_P2 = 2.0 * y[0] * squared_sum * (squared_sum + 1.0) / pow(1.0 + squared_sum * squared_sum + 2.0 * squared_dif, 1.5);
	double D4_v_P1 = 2.0 * y[1] / sqrt(1.0 + squared_sum * squared_sum + 2.0 * squared_dif);
	double D4_v_P2 = 2.0 * y[1] * squared_sum * (squared_sum - 1.0) / pow(1.0 + squared_sum * squared_sum + 2.0 * squared_dif, 1.5);
	double D4_u = D4_u_P1 - D4_u_P2;
	double D4_v = D4_v_P1 - D4_v_P2;
	double potential_u = 2.0 * (D1_u + 2.0 * (1.0 - mu) * D2_u + (1.0 - mu - C) * D3_u + 2.0 * (1.0 - mu) * D4_u);
	double potential_v = 2.0 * (D1_v + 2.0 * (1.0 - mu) * D2_v + (1.0 - mu - C) * D3_v + 2.0 * (1.0 - mu) * D4_v);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = potential_u + 8.0 * squared_sum * y[3];
	f[3] = potential_v - 8.0 * squared_sum * y[2];
	return GSL_SUCCESS;
}

int field_rotational(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = omega_x + 2.0 * y[3];
	f[3] = omega_y - 2.0 * y[2];
	return GSL_SUCCESS;
}

int jacobian_rotational(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double K1 = (1. - mu) / (r_1 * r_1 * r_1);
	double K2 = mu / (r_2 * r_2 * r_2);
	double A = (y[0] + mu) * (y[0] + mu) / (r_1 * r_1);
	double B = (y[0] + mu - 1.) * (y[0] + mu - 1.) / (r_2 * r_2);
	double C = y[1] * y[1] / (r_1 * r_1);
	double D = y[1] * y[1] / (r_2 * r_2);
	double E = (1. - mu) * (y[0] + mu) / (r_1 * r_1 * r_1 * r_1 * r_1);
	double F = mu * (y[0] + mu - 1.) / (r_2 * r_2 * r_2 * r_2 * r_2);
	double omega_xx = 1. - K1 * (1. - 3. * A) - K2 * (1. - 3. * B);
	double omega_yy = 1. - K1 * (1. - 3. * C) - K2 * (1. - 3. * D);
	double omega_xy = 3. * y[1] * (E + F);
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix *m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 0.0);
	gsl_matrix_set(m, 0, 2, 1.0);
	gsl_matrix_set(m, 0, 3, 0.0);
	gsl_matrix_set(m, 1, 0, 0.0);
	gsl_matrix_set(m, 1, 1, 0.0);
	gsl_matrix_set(m, 1, 2, 0.0);
	gsl_matrix_set(m, 1, 3, 1.0);
	gsl_matrix_set(m, 2, 0, omega_xx);
	gsl_matrix_set(m, 2, 1, omega_xy);
	gsl_matrix_set(m, 2, 2, 0.0);
	gsl_matrix_set(m, 2, 3, 2.0);
	gsl_matrix_set(m, 3, 0, omega_xy);
	gsl_matrix_set(m, 3, 1, omega_yy);
	gsl_matrix_set(m, 3, 2, -2.0);
	gsl_matrix_set(m, 3, 3, 0.0);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;
	return GSL_SUCCESS;
}