set terminal png
set ylabel "speedup"
set xlabel "number of cores"
set xtics 1
set title "nbody_brute_force_async Speedup with OpenMP execution"
set output "./MPI_nbody_brute_force_async_22_10_22_1666400736/SpeedUp-MPI.png"
plot x title 'Speedup max' with lines, './MPI_nbody_brute_force_async_22_10_22_1666400736/res_22_10_22_1666400736.data' using 1:($2/$3) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines, './MPI_nbody_brute_force_async_22_10_22_1666400736/res_22_10_22_1666400736.data' using 1:($2/$3/($1)) title 'Efficacité' with linespoints;
