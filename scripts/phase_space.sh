reset

set terminal qt 1 size 600, 600 font "Helvetica,15"

set loadpath "../results/phase_space"

set xlabel "Intersection number"
set ylabel "Precision"
unset key
set title "Jacobi constant precision on the phase space x-~x{1.1.}"

plot 'phase_space_jacobi_constant_precision.dat' w d

reset

set terminal qt 2 size 600, 600 font "Helvetica,15"

set loadpath "../results/phase_space"

set xlabel "x"
set ylabel "~x{1.1.}" offset 0,1,0
set xrange[-0.85:1.2]
set yrange[-2.5:2.5]
unset key
set title "Phase space x-~x{1.1.}"

plot 'phase_space.dat' w d, \
     'phase_space_initial_conditions.dat' w p, \

pause -1