#ifndef AUX_FUNC_VMO_H
#define AUX_FUNC_VMO_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PBSTR "=================================================="
#define PBWIDTH 50

/**
 * personal library created by VMO to assist
 * writting c programming codes
**/

/*************************** makes life easier ********************************/

// defines m as a square matrix of size dim by dim
// with elelents copied from the array x of size dim squared
int array_to_matrix(double **m, double *x, int dim);

// defines x as a copy of y
// where x and y are arrays of dimension dim
int copy(double *x, double *y, int dim);

// defines an interval of x as a copy of the same interval of y
// the copied indices go from initial_value to final_value
int copy_parts(double *x, double *y, int initial_value, int final_value);

// prints array of doubles of size dim
int print_array(double *x, int dim);

// prints the progress of a calculation as a filling bar
void print_prog(double percentage);

// returns one if n is positive and minus one if n is negative
double sign(double x);

// divides the dim-dimensional line connecting
// x_initial to x_final in number pieces
// and storages them in x
int split_section(double **x, double *x_initial, double *x_final, int number, int dim);

// switches the i-th and j-th elements of an array x
int switch_array_element(double *x, int i, int j);

/****************************** mathematics ***********************************/

// return the norm of an real array x of dimension dim
double array_norm(double *x, int dim);

// returns the euclidean distance between
// arrays x and y of dimension dim
double dist(double *x, double *y, int dim);

// returns the euclidean distance between
// a range of indicies of arrays x and y
double dist_partial(double *x, double *y, int index_1, int index_2);

// defines x as a two-dimensional real-valued identity matrix
// of size order by order in array form starting from index
// arr_first_index
int identity_matrix_array_form(double *x, int arr_first_index, int order);

// defines x as a linear combination of y and z
// where x, y and z are arrays of dimension dim
// and lamb is a constant parameter
int lin_comb(double *x, double lamb, double *y, double *z, int dim);

// defines x as a linear combination of y and z
// where x, y, z and lamb are arrays of dimension dim
int lin_comb_vec(double *x, double *lamb, double *y, double *z, int dim);

// transposes a real two-dimensional square matrix
// of size dim by dim
int transpose_2d_square_matrix(double **m, int dim);

/**************************** memory handling *********************************/

// allocates memory for one-dimensional double array x of size n
int alloc_1d_double(double **x, int n);

// allocates memory for one-dimensional integer array x of size n
int alloc_1d_int(int **x, int n);

// allocates memory for two-dimensional double array x of size n by n
int alloc_2d_double(double ***x, int n, int m);

// allocates memory for two-dimensional integer array x of size n by n
int alloc_2d_int(int ***x, int n, int m);

// allocates memory for three-dimensional double array x of size n by n by n
int alloc_3d_double(double ****x, int n, int m, int p);

// deallocates memory for one-dimensional double array x of size n
int dealloc_1d_double(double **x);

// deallocates memory for one-dimensional integer array x of size n
int dealloc_1d_int(int **x);

// deallocates memory for two-dimensional double array x of size n by n
int dealloc_2d_double(double ***x, int n);

// deallocates memory for two-dimensional integer array x of size n by n
int dealloc_2d_int(int ***x, int n);

// deallocates memory for three-dimensional double array x of size n by n by n
int dealloc_3d_double(double ****x, int n, int m);

// reallocates memory for one-dimensional double array x of size n
int realloc_1d_double(double **x, int n);

// reallocates memory for two-dimensional double array x of size n by n
int realloc_2d_double(double ***x, int n, int m, int o);

#endif