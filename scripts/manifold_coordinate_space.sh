reset

set terminal qt 1 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/coordinate_space"

set xlabel "Time"
set ylabel "Precision"
unset key
set title "Jacobi constant precision on the\ninvariant manifolds of the Lyapunov orbit around L1"

plot 'jacobi_constant_precision_manifold_unstable_left_L1.dat' lc rgb "red" w d, \
     'jacobi_constant_precision_manifold_unstable_right_L1.dat' lc rgb "blue" w d, \

reset

set terminal qt 2 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/coordinate_space"

set xlabel "Time"
set ylabel "Precision"
unset key
set title "Jacobi constant precision on the\ninvariant manifolds of the Lyapunov orbit around L2"

plot 'jacobi_constant_precision_manifold_unstable_right_L2.dat' lc rgb "web-blue" w d, \
     'jacobi_constant_precision_manifold_unstable_left_L2.dat' lc rgb "orange-red" w d, \

reset

set terminal qt 3 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/coordinate_space"

set xlabel "x"
set ylabel "y"
unset key
set title "Invariant manifolds of the Lyapunov orbits around L1 and L2"

plot 'manifold_unstable_left_L2.dat' lc rgb "orange-red" w l, \
     'manifold_stable_left_L2.dat' lc rgb "web-blue" w l, \
     'manifold_unstable_right_L2.dat' lc rgb "orange-red" w l, \
     'manifold_stable_right_L2.dat' lc rgb "web-blue" w l, \
     'manifold_unstable_left_L1.dat' lc rgb "red" w l, \
     'manifold_stable_left_L1.dat' lc rgb "blue" w l, \
     'manifold_unstable_right_L1.dat' lc rgb "red" w l, \
     'manifold_stable_right_L1.dat' lc rgb "blue" w l, \
     '../../lyapunov_orbit/lyapunov_orbit_L1.dat' lc rgb "black" w l, \
     '../../lyapunov_orbit/lyapunov_orbit_L2.dat' lc rgb "black" w l, \
     '../../zvc/zvc.dat' lc rgb "gray20" w l

pause -1