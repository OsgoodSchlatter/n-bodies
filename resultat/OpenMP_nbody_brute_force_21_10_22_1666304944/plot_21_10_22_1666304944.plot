set terminal png
set ylabel "speedup"
set xlabel "number of cores"
set xtics 1
set title "nbody_brute_force Speedup with OpenMP execution"
set output "./OpenMP_nbody_brute_force_21_10_22_1666304944/SpeedUp-OpenMP.png"
plot x title 'Speedup max' with lines, './OpenMP_nbody_brute_force_21_10_22_1666304944/res_21_10_22_1666304944.data' using 1:($2/$3) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines, './OpenMP_nbody_brute_force_21_10_22_1666304944/res_21_10_22_1666304944.data' using 1:($2/$3/($1)) title 'Efficacité' with linespoints;
