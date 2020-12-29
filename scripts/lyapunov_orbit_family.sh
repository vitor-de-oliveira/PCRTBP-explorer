reset

set terminal qt 1 size 600, 600 font "Helvetica,15"

set loadpath "../results/lyapunov_orbit"

set xlabel "x"
set ylabel "y"
set xrange[0.7:1.3]
set yrange[-0.35:0.35]
unset key
set title "Family of Lyapunov orbits"

plot 'family_L1.dat' w l, \
     'family_L2.dat' w l, \
     '../location/masses_position.dat' lc rgb "brown" pt 7 ps 1.5 w p, \
     '../location/lagrangian_points_position.dat' lc rgb "black" pt 7 ps 1 w p, \

pause -1