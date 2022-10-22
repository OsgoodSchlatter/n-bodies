#!/bin/bash

algo=nbody_barnes_hut

N_PROCESS=4

N_PARTICULE=1000
T_FINAL=6

date=$(date +"%d_%m_%y_%s")
hostsfile="./hosts"
dirname=MPI_$algo\_$date

seq_duration=0

#Build the structure of data
#./$date
#./$date/log
#./$date/res_$date.data

mkdir $dirname
mkdir $dirname/log

#mpirun de chauffe
mpirun -n N_PROCESS ../MPI/$algo $N_PARTICULE $T_FINAL

#Param init
echo \#N_PROCESS N_PARTICULE T_FINAL > ./$dirname/param.data
echo $N_PROCESS $N_PARTICULE $T_FINAL > ./$dirname/param.data >&1

echo \#n_process t_seq t_parallel > ./$dirname/res_$date.data

../sequential/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_seq.log >&1
seq_duration=$(cat ./$dirname/log/log_seq.log | grep "Simulation" | cut -d " " -f 3)
echo  1 $seq_duration $seq_duration >> ./$dirname/res_$date.data

mpirun -n $N_PROCESS ../MPI/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_$N_PROCESS.log >&1
duration=$(cat ./$dirname/log/log_$N_PROCESS.log | grep "Simulation" | cut -d " " -f 3)
echo  $N_PROCESS $seq_duration $duration >> ./$dirname/res_$date.data

#AFFICHAGE PLOT
#$cp ./plot.plot ./$date/plot_$date.plot
echo set terminal png >> ./$dirname/plot_$date.plot
echo set ylabel \"speedup\" >> ./$dirname/plot_$date.plot
echo set xlabel \"number of cores\" >> ./$dirname/plot_$date.plot
echo set xtics 1 >> ./$dirname/plot_$date.plot
echo set title \"$algo Speedup with OpenMP execution\" >> ./$dirname/plot_$date.plot

echo set output \"./$dirname/SpeedUp-MPI.png\" >> ./$dirname/plot_$date.plot
echo "plot x title 'Speedup max' with lines," \'./$dirname/res_$date.data\' "using 1:(\$2/\$3) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines," \'./$dirname/res_$date.data\' "using 1:(\$2/\$3/(\$1)) title 'Efficacité' with linespoints;" >> ./$dirname/plot_$date.plot

gnuplot ./$dirname/plot_$date.plot
