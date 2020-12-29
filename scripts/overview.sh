reset

set terminal qt 1 size 600, 600 font "Helvetica,15"

set loadpath "../results/"

set xlabel "Time"
set ylabel "Precision"
unset key
set title "Precision of the Jacobi constant on the traced orbit"

plot 'orbit/orbit_jacobi_constant_precision.dat' w d, \

reset

set terminal qt 2 size 600, 600 font "Helvetica,15"

set loadpath "../results/"

set xlabel "x"
set ylabel "y"
unset key
set title "Overview of the system"

plot 'orbit/orbit.dat' w l, \
     'location/masses_position.dat' lc rgb "brown" pt 7 ps 1.5 w p, \
     'location/lagrangian_points_position.dat' lc rgb "black" pt 7 ps 1 w p, \
     'zvc/zvc.dat' lc rgb "gray20" w l

pause -1