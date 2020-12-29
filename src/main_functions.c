#include "main_functions.h"

int trace_orbit(double *y, void *params, double tf)
{
	// declare and open exit files
	FILE *orb = fopen("results/orbit/orbit.dat", "w");
	FILE *ene = fopen("results/orbit/orbit_jacobi_constant_precision.dat", "w");

	// declare variables
	int step_number, size;
	double **orbit, *time;
	double J, J0;

	J0 = motion_integral(y, params, 0.0, "rotational");

	// evolve system
	evolve(y, params, tf, &orbit, &time, &step_number, &size);

	// write orbit and constant error to file
	for (int i = 0; i < step_number; i++)
	{
		fprintf(orb, "%1.15e %1.15e\n", orbit[i][0], orbit[i][1]);
		J = motion_integral(orbit[i], params, 0.0, "rotational");
		fprintf(ene, "%1.15e %1.15e\n", time[i], fabs(J-J0));
	}

	// close files
	fclose(orb);
	fclose(ene);

	// free memory
	dealloc_1d_double(&time);
	dealloc_2d_double(&orbit, size);

	printf("Data written in folder results/orbit\n");

	return 0;
}

int draw_phase_space(double motion_constant, void *params, double tf, double coordinate_min, double coordinate_max, double velocity_min, double velocity_max, int nc, int nv)
{
	// declare and open exit files
	FILE *psp = fopen("results/phase_space/phase_space.dat", "w");
	FILE *inc = fopen("results/phase_space/phase_space_initial_conditions.dat", "w");
	FILE *ene = fopen("results/phase_space/phase_space_jacobi_constant_precision.dat", "w");

	// declare variables
	int cross_coordinate, cross_velocity;
	int map_coordinate, map_velocity;
	int size_orbit, size_map;
	int step_number, step_number_map;
	int size_orbit_bw, size_map_bw;
	int step_number_bw, step_number_map_bw;
	double y[4], y0[4];
	double coordinate, velocity, cross_condition;
	double **orbit, *time, **map, *energy;
	double **orbit_bw, *time_bw, **map_bw, *energy_bw;

	// initialize ic parameters
	// same as map definition
	cross_condition = 0.0;
	cross_coordinate = 1;
	cross_velocity = 3;
	map_coordinate = 0;
	map_velocity = 2;

	// loop over coordinate values
	coordinate = coordinate_min;
	for (int i = 0; i < nc; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", i + 1, nc);

		// loop over velocity values
		velocity = velocity_min;
		for (int j = 0; j < nv; j++)
		{
			// print progress on velocity
			print_prog((double)(j + 1) / (double)nv);

			// calculate initial condition
			y[map_coordinate] = coordinate;
			y[cross_coordinate] = cross_condition;
			y[map_velocity] = velocity;
			integral_to_velocity(y, params, 0.0, motion_constant, "rotational", cross_velocity);

			// check if valid initial condition
			if (y[cross_velocity] != y[cross_velocity])
			{
				// printf("Warning: initial value outside range\n");
				goto final;
			}

			// write initial condition to file
			fprintf(inc, "%1.15e %1.15e\n", coordinate, velocity);

			// store initial condition
			copy(y0, y, 4);

			// calculate forward orbit
			evolve(y, params, tf, &orbit, &time, &step_number, &size_orbit);

			// calculate forward map
			poincare_map(params, &map, &energy, NULL, &step_number_map, &size_map, orbit, time, step_number);

			// write forward map to file
			for (int k = 0; k < step_number_map; k++)
			{
				fprintf(psp, "%1.15e %1.15e\n", map[k][0], map[k][1]);
				fprintf(ene, "%d %1.15e\n", k, fabs(energy[k] - motion_constant));
			}
			
			// free memory forward
			dealloc_1d_double(&time);
			dealloc_2d_double(&orbit, size_orbit);
			dealloc_1d_double(&energy);
			dealloc_2d_double(&map, size_map);
			
			// calculate backward orbit
			evolve(y0, params, -tf, &orbit_bw, &time_bw, &step_number_bw, &size_orbit_bw);

			// calculate backward map
			poincare_map(params, &map_bw, &energy_bw, NULL, &step_number_map_bw, &size_map_bw, orbit_bw, time_bw, step_number_bw);
			
			// write backward map to file
			for (int k = 0; k < step_number_map_bw; k++)
			{
				fprintf(psp, "%1.15e %1.15e\n", map_bw[k][0], map_bw[k][1]);
				fprintf(ene, "%d %1.15e\n", k, fabs(energy_bw[k] - motion_constant));
			}

			// free memory backward
			dealloc_1d_double(&time_bw);
			dealloc_2d_double(&orbit_bw, size_orbit_bw);
			dealloc_1d_double(&energy_bw);
			dealloc_2d_double(&map_bw, size_map_bw);

			// create new line on exit file
			fprintf(psp, "\n");
			fprintf(ene, "\n");

		// jump here if invalid initial condition
		final:;
			
			// update valocity
			if (nv > 1)
			{
				velocity += fabs(velocity_max - velocity_min) / (double)(nv - 1);
			}
		}

		// update coordinate
		if (nc > 1)
		{
			coordinate += fabs(coordinate_max - coordinate_min) / (double)(nc - 1);
		}
		
		// new line on terminal
		printf("\n");
	}

	// close exit files
	fclose(psp);
	fclose(inc);
	fclose(ene);

	printf("Data written in folder results/phase_space\n");

	return 0;
}

int trace_Lyapunov_orbit(void *params, double J, int n)
{
	// declare variables
	double y[4], period;

	// find orbit period and initial condition
	find_Lyapunov_orbit_parameters(params, J, y, &period, n);

	// declare and open exit files
	char filename1[100], filename2[100];
	sprintf(filename1, "results//lyapunov_orbit/lyapunov_orbit_L%d.dat", n);
	sprintf(filename2, "results//lyapunov_orbit/jacobi_constant_lyapunov_orbit_L%d.dat", n);
	FILE *orb = fopen(filename1, "w");
	FILE *ene = fopen(filename2, "w");

	// declare variables
	int step_number, size;
	double **orbit, *time;

	// evolve system
	evolve(y, params, period, &orbit, &time, &step_number, &size);

	// write orbit and constant error to file
	for (int i = 0; i < step_number; i++)
	{
		fprintf(orb, "%1.15e %1.15e\n", orbit[i][0], orbit[i][1]);
		J = motion_integral(orbit[i], params, 0.0, "rotational");
		fprintf(ene, "%1.15e %1.15e\n", time[i], J);
	}

	// close files
	fclose(orb);
	fclose(ene);

	// free memory
	dealloc_1d_double(&time);
	dealloc_2d_double(&orbit, size);

	printf("Data written in folder results/lyapunov_orbit\n");

	return 0;
}

int trace_Lyapunov_orbit_family(void *params, int number_of_orbits, double delta_J, int n)
{
	// declare and open exit files
	char filename[100];
	sprintf(filename, "results/lyapunov_orbit/family_L%d.dat", n);
	FILE *fam = fopen(filename, "w");

	// declare variables
	int step_number, size;
	double y[4], period, J;
	double **orbit, *time;
	double constant_precision;
	double y0[4], y01[4], y02[4], yk[4], y03;
	double half_period_try, Jupo, periodupo;
	double delta_x, step_x, periodk;
	double delta_ydot, delta_ydotk;
	double ykeep[4], delta_ydot_keep, periodkeep;
	bool keep_values;

	// define method parameters
	constant_precision = 1e-13;
	delta_x = 1e-4;
	step_x = 10. * delta_x;

	// calculate reference orbits in linear region
	Lyapunov_orbit_initial_condition_linear_region(params, y01, delta_x, &Jupo, &periodupo, n);
	Lyapunov_orbit_initial_condition_linear_region(params, y02, step_x, &Jupo, &periodupo, n);

	// first parameters
	delta_ydot = fabs(y02[3] - y01[3]);
	copy(y0, y02, 4);

	// jacobi constant of n-th Lagrangian point
	J = lagrangian_point_jacobi_constant(params, n);

	// loop for family members
	for (int i = 0; i < number_of_orbits; i++)
	{
		// jacobi for family member
		J -= delta_J;

		printf("Jacobi constant = %1.5e\n", J);

		// start method
		keep_values = true;
		do
		{
			do
			{
				y03 = y0[3];
				copy(yk, y0, 4);
				delta_ydotk = delta_ydot;
				periodk = periodupo;
				y0[0] += step_x;
				y0[1] = 0.0;
				y0[2] = 0.0;
				y0[3] -= delta_ydot;
				half_period_try = periodupo / 2.0;
				find_periodic_orbit_initial_parameters(params, y0, half_period_try, &Jupo, &periodupo);
				delta_ydot = fabs(y0[3] - y03);
			} while (J < Jupo);
			if (keep_values == true)
			{
				delta_ydot_keep = delta_ydot;
				copy(ykeep, yk, 4);
				periodkeep = periodk;
				keep_values = false;
			}
			copy(y0, yk, 4);
			periodupo = periodk;
			step_x /= 2.0;
			delta_ydot = delta_ydotk / 2.0;
			printf("Constant precision = %1.15e\n", fabs(J - Jupo));
		} while (fabs(J - Jupo) > constant_precision);

		// update initial condition and period of orbit
		copy(y, y0, 4);
		period = periodupo;
		
		// calculate member
		evolve(y, params, period, &orbit, &time, &step_number, &size);

		// write member to file
		for (int j = 0; j < step_number; j++)
		{
			fprintf(fam, "%1.15e %1.15e\n", orbit[j][0], orbit[j][1]);
		}
		fprintf(fam, "\n");

		// free memory
		dealloc_1d_double(&time);
		dealloc_2d_double(&orbit, size);

		// keep values for next member
		delta_ydot = delta_ydot_keep;
		copy(y0, ykeep, 4);
		periodupo = periodkeep;
		step_x = 10. * delta_x;
	}

	// close exit file
	fclose(fam);

	printf("Data written in folder results/lyapunov_orbit\n");

	return 0;
}

int trace_zvc(double J, void *params)
{
	FILE *orb = fopen("results/zvc/zvc.dat", "w");
	FILE *ene = fopen("results/zvc/zvc_energy.dat", "w");
	int status_root, status_zvc;
	double mu, J_step, t, tf, h, delta, x1, x2;
	double y0[2], y_zvc[2], y_J[4];
	mu = *(double *)params;
	delta = 1e-1;
	h = 1e-3;
	tf = 20.;
	if (J > lagrangian_point_jacobi_constant(params, 1))
	{
		printf("Warning: zero-velocity curve calculation\nwas not implemented yet\nfor this Jacobian constant");
		return 1;
	}
	else if (J >= lagrangian_point_jacobi_constant(params, 3))
	{
		printf("Motion constant between J1 and J3\n");
		x1 = collinear_lagrangian_point_x_coordinate(params, 3) + delta;
		x2 = collinear_lagrangian_point_x_coordinate(params, 3) - delta;
	}
	else if (J >= lagrangian_point_jacobi_constant(params, 4))
	{
		printf("Motion constant between J3 and J4\n");
		x1 = 0.5 * sqrt(3.0) + delta;
		x2 = -0.5 * sqrt(3.0) - delta;
	}
	else
	{
		printf("Warning: there is no zero-velocity curve\nfor this Jacobian constant\n");
		return 1;
	}
	const gsl_odeiv2_step_type *T_zvc = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_system sys_zvc = {field_zero_velociy_curves_rotational, NULL, 2, &mu};
	gsl_odeiv2_driver *d_zvc = gsl_odeiv2_driver_alloc_standard_new(&sys_zvc, T_zvc, h, 0.0, 1e-15, 1.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d_zvc, 1e-3);
	status_root = root_zvc(J, &mu, y0, x1);
	if (status_root != 0)
	{
		printf("Warning: could not trace zero-velocity curve\n");
		return 1;
	}
	copy(y_zvc, y0, 2);
	fprintf(orb, "%1.15e %1.15e\n", y_zvc[0], y_zvc[1]);
	t = 0.0;
	while (t < tf)
	{
		status_zvc = gsl_odeiv2_driver_apply (d_zvc, &t, t + sign(tf) * 1e-3, y_zvc);
		if (status_zvc != GSL_SUCCESS)
		{
			printf("Warning: error ZVC, %s\n", gsl_strerror(status_zvc));
			exit(2);
		}
		y_J[0] = y_zvc[0];
		y_J[1] = y_zvc[1];
		y_J[2] = 0.;
		y_J[3] = 0.;
		J_step = motion_integral(y_J, &mu, 0.0, "rotational");
		fprintf(ene, "%1.15e %1.15e\n", t, J_step);
		fprintf(orb, "%1.15e %1.15e\n", y_zvc[0], y_zvc[1]);
		if (t > 1. && (dist(y_zvc, y0, 2) < 1e-3))
		{
			break;
		}
	}
	status_root = root_zvc(J, &mu, y0, x2);
	if (status_root != 0)
	{
		printf("Warning: could not trace zero-velocity curve\n");
		return 1;
	}
	copy(y_zvc, y0, 2);
	fprintf(orb, "\n%1.15e %1.15e\n", y_zvc[0], y_zvc[1]);
	t = 0.0;
	while (t < tf)
	{
		status_zvc = gsl_odeiv2_driver_apply (d_zvc, &t, t + sign(tf) * 1e-3, y_zvc);
		if (status_zvc != GSL_SUCCESS)
		{
			printf("Warning: error ZVC, %s\n", gsl_strerror(status_zvc));
			exit(2);
		}
		y_J[0] = y_zvc[0];
		y_J[1] = y_zvc[1];
		y_J[2] = 0.;
		y_J[3] = 0.;
		J_step = motion_integral(y_J, &mu, 0.0, "rotational");
		fprintf(ene, "\n%1.15e %1.15e", t, J_step);
		fprintf(orb, "%1.15e %1.15e\n", y_zvc[0], y_zvc[1]);
		if (t > 1. && (dist(y_zvc, y0, 2) < 1e-3))
		{
			break;
		}
	}
	fclose(orb);
	fclose(ene);
	gsl_odeiv2_driver_free(d_zvc);

	printf("Data written in folder results/zvc\n");

	return 0;
}

int trace_manifolds_Lyapunov_orbit(void *params, double J, double time_left, double time_right, int man_ic_number, int n)
{
	// check if Lagrangian number is valid
	if (n != 1 && n != 2)
	{
		printf("Warning: invalid index\nfor Lyapunov orbit in\nfunction trace_manifolds_Lyapunov_orbit\n");
		printf("Warning: function not yet implemented\nfor L3, L4 and L5\n");
		exit(2);
	}

	// declare exit files
	FILE *mup, *mum, *msp, *msm;
	FILE *mupj, *mumj;
	FILE *mupm, *mumm, *mspm, *msmm;
	FILE *mupmj, *mummj;
	char filename[100];

	// open exit files

	// manifolds on real space
	sprintf(filename, "results/manifolds/coordinate_space/manifold_unstable_left_L%d.dat", n);
	mup = fopen(filename, "w");

	sprintf(filename, "results/manifolds/coordinate_space/manifold_unstable_right_L%d.dat", n);
	mum = fopen(filename, "w");

	sprintf(filename, "results/manifolds/coordinate_space/manifold_stable_left_L%d.dat", n);
	msp = fopen(filename, "w");

	sprintf(filename, "results/manifolds/coordinate_space/manifold_stable_right_L%d.dat", n);
	msm = fopen(filename, "w");

	// jacobi constant of 
	// manifolds on real space
	sprintf(filename, "results/manifolds/coordinate_space/jacobi_constant_precision_manifold_unstable_left_L%d.dat", n);
	mupj = fopen(filename, "w");

	sprintf(filename, "results/manifolds/coordinate_space/jacobi_constant_precision_manifold_unstable_right_L%d.dat", n);
	mumj = fopen(filename, "w");

	// manifolds on poincare map
	sprintf(filename, "results/manifolds/phase_space/manifold_unstable_left_map_L%d.dat", n);
	mupm = fopen(filename, "w");

	sprintf(filename, "results/manifolds/phase_space/manifold_unstable_right_map_L%d.dat", n);
	mumm = fopen(filename, "w");

	sprintf(filename, "results/manifolds/phase_space/manifold_stable_left_map_L%d.dat", n);
	mspm = fopen(filename, "w");

	sprintf(filename, "results/manifolds/phase_space/manifold_stable_right_map_L%d.dat", n);	
	msmm = fopen(filename, "w");

	// jacobi constant of 
	// manifolds on poincare map
	sprintf(filename, "results/manifolds/phase_space/jacobi_constant_precision_manifold_unstable_left_map_L%d.dat", n);	
	mupmj = fopen(filename, "w");

	sprintf(filename, "results/manifolds/phase_space/jacobi_constant_precision_manifold_unstable_right_map_L%d.dat", n);
	mummj = fopen(filename, "w");

	// declare variables
	int count;
	int size_orbit_mum, size_orbit_mup;
	int size_map_mumm, size_map_mupm;
	int step_number_mum, step_number_mup;
	int step_number_mumm, step_number_mupm;
	int *map_size, *new_map_size;	
	int *crossing_min;
	double y0[4], y[20], **man_ic, ***map;
	double man_ic_dist_from_upo;
	double period, man_ic_step, norm, J_man;
	double monodromy_matrix[16], eig_val[2], eig_vec[2][4];
	double eig_val_real[4], eig_val_imag[4];
	double eig_vec_real[4][4], eig_vec_imag[4][4];
	double transition_matrix[16], eig_vec_uns_evolve[4];
	double man_ic_uns_plus[4], man_ic_uns_minus[4];
	double **orbit_mum, *time_mum, **map_mum, *map_energy_mum;
	double **orbit_mup, *time_mup, **map_mup, *map_energy_mup;
	double t_man_plus, t_man_minus;

	// define function parameters
	man_ic_dist_from_upo = 1e-6;

	// allocate memory
	alloc_1d_int(&map_size, man_ic_number);
	alloc_1d_int(&new_map_size, man_ic_number);
	alloc_1d_int(&crossing_min, man_ic_number);
	alloc_2d_double(&man_ic, man_ic_number, 20);
	alloc_3d_double(&map, man_ic_number, 50, 2);

	// calculate Lyapunov orbit
	find_Lyapunov_orbit_parameters(params, J, y0, &period, n);
	copy(y, y0, 4);
	identity_matrix_array_form(y, 4, 4);

	// determine distance between 
	// points of Lyapunov orbit
	man_ic_step = period / (double)man_ic_number;

	// calculate points on Lyapunov orbit
	// along the transition matrices
	for (int i = 0; i < man_ic_number; i++)
	{
		evolve_extended(y, params, man_ic_step);
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
		exit(2);
	}

	// defining branches
	if (eig_vec[0][0] > 0.0)
	{
		t_man_plus = time_right;
		t_man_minus = time_left;
	}
	else
	{
		t_man_plus = time_left;
		t_man_minus = time_right;
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
		evolve(man_ic_uns_minus, params, t_man_minus, &orbit_mum, &time_mum, &step_number_mum, &size_orbit_mum);

		// write manifold orbit minus branch
		for (int j = 0; j < step_number_mum; j++)
		{
			// write unstable manifold
			fprintf(mum, "%1.15e %1.15e %1.15e\n", orbit_mum[j][0], orbit_mum[j][1], orbit_mum[j][2]);

			// write stable manifold
			fprintf(msm, "%1.15e %1.15e %1.15e\n", orbit_mum[j][0], -orbit_mum[j][1], -orbit_mum[j][2]);

			// calculate jacobi constant
			J_man = motion_integral(orbit_mum[j], params, 0.0, "rotational");

			// write manifold orbit error
			fprintf(mumj, "%1.15e %1.15e\n", time_mum[j], fabs(J_man - J));

		}
		fprintf(mum, "\n");
		fprintf(msm, "\n");
		fprintf(mumj, "\n");

		// calculate manifold orbit on poincare map
		poincare_map(params, &map_mum, &map_energy_mum, NULL, &step_number_mumm, &size_map_mumm, orbit_mum, time_mum, step_number_mum);
		
		// keep orbit sizes
		map_size[i] = step_number_mumm;

		// write manifold orbit minus branch on poincare map
		for (int j = 0; j < step_number_mumm; j++)
		{
			// keep map values
			copy(map[i][j], map_mum[j], 2);
			
			// write unstable manifold on poincare map
			fprintf(mumm, "%1.15e %1.15e\n", map_mum[j][0], map_mum[j][1]);

			// write stable manifold on poincare map
			fprintf(msmm, "%1.15e %1.15e\n", map_mum[j][0], -map_mum[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mummj, "%d %1.15e\n", j, fabs(map_energy_mum[j] - J_man));

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
		evolve(man_ic_uns_plus, params, t_man_plus, &orbit_mup, &time_mup, &step_number_mup, &size_orbit_mup);

		// write manifold orbit plus branch
		for (int j = 0; j < step_number_mup; j++)
		{
			// write unstable manifold
			fprintf(mup, "%1.15e %1.15e %1.15e\n", orbit_mup[j][0], orbit_mup[j][1], orbit_mup[j][2]);

			// write stable manifold
			fprintf(msp, "%1.15e %1.15e %1.15e\n", orbit_mup[j][0], -orbit_mup[j][1], -orbit_mup[j][2]);

			// calculate jacobi constant
			J_man = motion_integral(orbit_mup[j], params, 0.0, "rotational");
			
			// write manifold orbit error
			fprintf(mupj, "%1.15e %1.15e\n", time_mup[j], fabs(J_man - J));
		}
		fprintf(mup, "\n");
		fprintf(msp, "\n");
		fprintf(mupj, "\n");

		// calculate manifold orbit on poincare map
		poincare_map(params, &map_mup, &map_energy_mup, NULL, &step_number_mupm, &size_map_mupm, orbit_mup, time_mup, step_number_mup);

		// write manifold orbit plus branch on poincare map
		for (int j = 0; j < step_number_mupm; j++)
		{
			// write unstable manifold on poincare map
			fprintf(mupm, "%1.15e %1.15e\n", map_mup[j][0], map_mup[j][1]);

			// write stable manifold on poincare map
			fprintf(mspm, "%1.15e %1.15e\n", map_mup[j][0], -map_mup[j][1]);

			// write manifold orbit error on poincare map
			fprintf(mupmj, "%d %1.15e\n", j, fabs(map_energy_mup[j] - J_man));
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
	dealloc_1d_int(&map_size);
	dealloc_1d_int(&new_map_size);
	dealloc_1d_int(&crossing_min);
	dealloc_2d_double(&man_ic, man_ic_number);
	dealloc_3d_double(&map, man_ic_number, 50);

	// close exit files
	fclose(mup); fclose(mum); fclose(msp); fclose(msm);
	fclose(mupj); fclose(mumj);
	fclose(mupm); fclose(mumm); fclose(mspm); fclose(msmm);
	fclose(mupmj); fclose(mummj);

	printf("Data written in folder results/manifolds/space\n");
	printf("Data written in folder results/manifolds/map\n");

	return 0;
}
