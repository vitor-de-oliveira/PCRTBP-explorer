reset

set terminal qt 1 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/phase_space"

set xlabel "Intersection number"
set ylabel "Precision"
unset key
set title "Jacobi constant precision on the invariant manifolds \n of the Lyapunov orbit around L1 measured \n in the phase space"

plot 'jacobi_constant_precision_manifold_unstable_left_map_L1.dat' lc rgb "red" pt 7 ps 0.5 w p, \
     'jacobi_constant_precision_manifold_unstable_right_map_L1.dat' lc rgb "blue" pt 7 ps 0.5 w p, \

reset

set terminal qt 2 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/phase_space"

set xlabel "Intersection number"
set ylabel "Precision"
unset key
set title "Jacobi constant precision on the invariant manifolds \n of the Lyapunov orbit around L2 measured \n in the phase space"

plot 'jacobi_constant_precision_manifold_unstable_left_map_L2.dat' lc rgb "orange-red" pt 7 ps 0.5 w p, \
     'jacobi_constant_precision_manifold_unstable_right_map_L2.dat' lc rgb "web-blue" pt 7 ps 0.5 w p, \

reset

set terminal qt 3 size 600, 600 font "Helvetica,15"

set loadpath "../results/manifolds/phase_space"

set xlabel "x"
set ylabel "~x{1.1.}"
unset key
set title "Invariant manifolds of the Lyapunov orbits \n around L1 and L2 in the phase space x-~x{1.1.}"

plot 'manifold_unstable_left_map_L2.dat' lc rgb "orange-red" pt 7 ps 0.5 w p, \
     'manifold_stable_left_map_L2.dat' lc rgb "web-blue" pt 7 ps 0.5 w p, \
     'manifold_unstable_right_map_L2.dat' lc rgb "orange-red" pt 7 ps 0.5 w p, \
     'manifold_stable_right_map_L2.dat' lc rgb "web-blue" pt 7 ps 0.5 w p, \
     'manifold_unstable_left_map_L1.dat' lc rgb "red" pt 7 ps 0.5 w p, \
     'manifold_stable_left_map_L1.dat' lc rgb "blue" pt 7 ps 0.5 w p, \
     'manifold_unstable_right_map_L1.dat' lc rgb "red" pt 7 ps 0.5 w p, \
     'manifold_stable_right_map_L1.dat' lc rgb "blue" pt 7 ps 0.5 w p, \

pause -1