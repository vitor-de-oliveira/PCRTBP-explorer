#include "auxiliar_functions_gsl.h"

int set_driver(gsl_odeiv2_driver *d, gsl_odeiv2_system *sys, gsl_odeiv2_step_type *T, char *error, char *control)
{
	if (strcmp(error, "absolute") == 0)
	{
		d = gsl_odeiv2_driver_alloc_standard_new(sys, T, 1e-3, 1e-15, 0.0, 0.0, 0.0);
	}
	else if (strcmp(error, "relative") == 0)
	{
		d = gsl_odeiv2_driver_alloc_standard_new(sys, T, 1e-3, 0.0, 1e-14, 1.0, 0.0);
	}
	else
	{
		printf("Warning: undefined error type\n");
		exit(2);
	}
	if (strcmp(control, "adaptive") == 0)
	{
		gsl_odeiv2_driver_set_hmax(d, 1e-3);
		gsl_odeiv2_driver_set_hmin(d, 1e-11);
	} 
	else if (strcmp(control, "fixed") != 0)
	{
		printf("Warning: undefined control type\n");
		exit(2);
	}
	return 0;
}

int set_integrator(const gsl_odeiv2_step_type *T, char *integrator)
{
    if (strcmp(integrator, "rk2") == 0)
	{
		T = gsl_odeiv2_step_rk2;
	}
	else if (strcmp(integrator, "rk4") == 0)
	{
		T = gsl_odeiv2_step_rk4;
	}
	else if (strcmp(integrator, "rk45") == 0)
	{
		T = gsl_odeiv2_step_rkf45;
	}
	else if (strcmp(integrator, "rkck") == 0)
	{
		T = gsl_odeiv2_step_rkck;
	}
	else if (strcmp(integrator, "rk8pd") == 0)
	{
		T = gsl_odeiv2_step_rk8pd;
	}
	else if (strcmp(integrator, "rk4imp") == 0)
	{
		T = gsl_odeiv2_step_rk4imp;
	}
	else if (strcmp(integrator, "bsimp") == 0)
	{
		T = gsl_odeiv2_step_bsimp;
	}
	else if (strcmp(integrator, "msadams") == 0)
	{
		T = gsl_odeiv2_step_msadams;
	}
	else if (strcmp(integrator, "msbdf") == 0)
	{
		T = gsl_odeiv2_step_msbdf;	
	}
	else
	{
		printf("Warning: invalid GSL integrator integrator\n");
		exit(2);
	}
	return 0;
}

int set_system(gsl_odeiv2_system *sys, void *par, char *system)
{
	double mu = *(double *)par;

	if (strcmp(system, "rotational") == 0 || 
		strcmp(system, "regularized") == 0)
	{
		(*sys).function = field_rotational;
		(*sys).jacobian = jacobian_rotational;
		(*sys).dimension = 4;
		(*sys).params = &mu;
	}
	else if (strcmp(system, "extended") == 0)
	{
		(*sys).function = field_extended;
		(*sys).jacobian = NULL;
		(*sys).dimension = 20;
		(*sys).params = &mu;
	}
	else if (strcmp(system, "hamiltonian") == 0)
	{
		(*sys).function = field_hamiltonian;
		(*sys).jacobian = NULL;
		(*sys).dimension = 4;
		(*sys).params = &mu;
	}
	else if (strcmp(system, "inertial") == 0)
	{
		(*sys).function = field_inertial;
		(*sys).jacobian = NULL;
		(*sys).dimension = 4;
		(*sys).params = &mu;
	}
	else
	{
		printf("Warning: undefined system\n");
		exit(2);
	}
	return 0;
}