set terminal png
set ylabel "speedup"
set xlabel "number of cores"
set xtics 1
set title "nbody_brute_force Speedup with MPI + OpenMP execution"
set output "./MPI_OpenMP_nbody_brute_force_21_10_22_1666308438/SpeedUp-MPI_OpenMP.png"
plot x title 'Speedup max' with lines, './MPI_OpenMP_nbody_brute_force_21_10_22_1666308438/res_21_10_22_1666308438.data' using 3:($4/$5) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines, './MPI_OpenMP_nbody_brute_force_21_10_22_1666308438/res_21_10_22_1666308438.data' using 3:($4/$5/($3)) title 'Efficacité' with linespoints;
