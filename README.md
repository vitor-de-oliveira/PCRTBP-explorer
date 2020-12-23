This software was developed to numerically investigate the dynamical and geometrical aspects of the planar circular restricted three-body problem. It is entirely written in c and most of the functions use the [Gnu Scientific Library](https://www.gnu.org/software/gsl/). In order to use it, open the main.c file and uncomment the function you want to use. You might have to input some values such as the Jacobi constant and the mass parameter. The results should appear as a .dat file in the results folder. The graphical part is up to you (for now).

The available tools include: numerical integration of solutions, using the methods provided by GSL; determination of periodic orbits, including Lyapunov and g-family orbits; calculation of two-dimensional invariant manifolds of unstable periodic orbits; tracing zero-velocity curves; determination of Lagrangian points location; and drawing of Poincare maps.

As a default: the system's reference frame is rotational, with the smaller mass to the right of the origin; the numerical integrator is the Prince-Dormand Runge-Kutta 8(9); the Poincare map is given by y=0 and ydot>0; the equations of motion are locally regularized around the primaries using Levi-Civita transformation; and collisions with the primaries are not considered. All of these settings can be changed in the source code.

As an example of the software capabilities, see the following paper: V. M. de Oliveira, P. A. Sousa-Silva, I. L. Caldas, "Order-chaos-order and invariant manifolds in the bounded planar Earth-Moon system," *Celestial Mechanics and Dynamical Astronomy*, vol. 132, no. 11, pp 1-17, 2020. DOI: 10.1007/s10569-020-09989-x.

Author: V. de Oliveira

Last update: Dec. 23 of 2020

## Requirements
```sh
* gcc
* GSL
* some knowledge of c
```

## Compile
```sh
make -s
```

## Run
```sh
./3BP
```

## Clean
```sh
make clean -s
```
