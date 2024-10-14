set terminal pngcairo enhanced size 800,900
set output 'plot.png'

set multiplot layout 3,1

set style line 1 lt rgb "red" lw 1
set style line 2 lt rgb "blue" lw 1

# Первый график: плотность
set ylabel "Density"
set grid
plot 'solution.dat' using 1:2 with lines ls 1 title "Godunov", \
    'exact_solution.dat' using 1:2 with lines ls 2 title "Exact"

# Второй график: скорость
set ylabel "Velocity"
set grid
plot 'solution.dat' using 1:3 with lines ls 1 title "Godunov", \
    'exact_solution.dat' using 1:3 with lines ls 2 title "Exact"

# Третий график: давление
set ylabel "Pressure"
set grid
plot 'solution.dat' using 1:4 with lines ls 1 title "Godunov", \
    'exact_solution.dat' using 1:4 with lines ls 2 title "Exact"
unset multiplot