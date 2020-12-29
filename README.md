This software was developed to numerically investigate the dynamical and geometrical aspects of the planar circular restricted three-body problem. It is entirely written in c and most of the functions use the [Gnu Scientific Library](https://www.gnu.org/software/gsl/). In order to use it, open the main.c file and uncomment the function you want. You might have to input some values such as the Jacobi constant and the mass parameter. The results should appear as a .dat file in the results folder. If gnuplot is installed, there are some shell scripts that can be used to plot the data files in the scripts folder.

The available tools include: numerical integration of solutions, using the methods provided by GSL; determination of periodic orbits, including Lyapunov and g-family orbits; calculation of two-dimensional invariant manifolds of unstable periodic orbits; tracing zero-velocity curves; determination of Lagrangian points location; and drawing of Poincare maps.

As a default: the system's reference frame is rotational, with the smaller mass to the right of the origin; the numerical integrator is the Prince-Dormand Runge-Kutta 8(9); the Poincare map is given by y=0 and ydot>0; the equations of motion are locally regularized around the primaries using Levi-Civita transformation; and collisions with the primaries are not considered. All of these settings can be changed in the source code.

As an example of the software capabilities, see the following paper: V. M. de Oliveira, P. A. Sousa-Silva, I. L. Caldas, "Order-chaos-order and invariant manifolds in the bounded planar Earth-Moon system," *Celestial Mechanics and Dynamical Astronomy*, vol. 132, no. 11, pp 1-17, 2020. DOI: 10.1007/s10569-020-09989-x.

Author: V. de Oliveira

Last update: Dec. 29 of 2020

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

## Shell scripts
```sh
gnuplot file.sh
```

![fig_overview](https://user-images.githubusercontent.com/52892492/103288957-d754ab00-49c4-11eb-9df5-57f29a20ed4b.png)

![fig_phase_space](https://user-images.githubusercontent.com/52892492/103288976-e0457c80-49c4-11eb-8cd2-1de7b88fb316.png)

![fig_manifolds_coordinate_space](https://user-images.githubusercontent.com/52892492/103288983-e3d90380-49c4-11eb-9562-15bdb4a4a7c3.png)

![fig_manifolds_coordinate_space_zoom](https://user-images.githubusercontent.com/52892492/103288986-e63b5d80-49c4-11eb-9258-3e972953ab68.png)
