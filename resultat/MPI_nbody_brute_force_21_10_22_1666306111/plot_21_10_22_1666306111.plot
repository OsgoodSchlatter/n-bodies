set terminal png
set ylabel "speedup"
set xlabel "number of cores"
set xtics 1
set title "nbody_brute_force Speedup with OpenMP execution"
set output "./MPI_nbody_brute_force_21_10_22_1666306111/SpeedUp-OpenMP.png"
plot x title 'Speedup max' with lines, './MPI_nbody_brute_force_21_10_22_1666306111/res_21_10_22_1666306111.data' using 1:($2/$3) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines, './MPI_nbody_brute_force_21_10_22_1666306111/res_21_10_22_1666306111.data' using 1:($2/$3/($1)) title 'Efficacité' with linespoints;
