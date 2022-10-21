#!/bin/bash

algo=nbody_brute_force

MAX_THREAD=48

N_PARTICULE=2000
T_FINAL=3

date=$(date +"%d_%m_%y_%s")
dirname=OpenMP_$algo\_$date

seq_duration=0

#Build the structure of data
#./$date
#./$date/log
#./$date/res_$date.data

mkdir $dirname
mkdir $dirname/log

#openmp de chauffe
OMP_NUM_THREADS=2 ../OpenMP/$algo $N_PARTICULE $T_FINAL

#Param init
echo \#MAX_THREAD N_PARTICULE T_FINAL > ./$dirname/param.data
echo $MAX_THREAD $N_PARTICULE $T_FINAL > ./$dirname/param.data >&1

echo \#n_thread t_seq t_parallel > ./$dirname/res_$date.data

#echo $n_process $MAX_THREAD
OMP_NUM_THREADS=1 ../OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_seq.log >&1

seq_duration=$(cat ./$dirname/log/log_seq.log | grep "Simulation" | cut -d " " -f 3)

echo  1 $seq_duration $seq_duration >> ./$dirname/res_$date.data

for n_thread in $(seq 2 $MAX_THREAD)
do
    OMP_NUM_THREADS=$n_thread ../OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_$n_thread.log >&1

    duration=$(cat ./$dirname/log/log_$n_thread.log | grep "Simulation" | cut -d " " -f 3)

    echo $n_thread $seq_duration $duration >> ./$dirname/res_$date.data
    echo $n_thread
done

#AFFICHAGE PLOT
#$cp ./plot.plot ./$date/plot_$date.plot
echo set terminal png >> ./$dirname/plot_$date.plot
echo set ylabel \"speedup\" >> ./$dirname/plot_$date.plot
echo set xlabel \"number of cores\" >> ./$dirname/plot_$date.plot
echo set xtics 1 >> ./$dirname/plot_$date.plot
echo set title \"$algo Speedup with OpenMP execution\" >> ./$dirname/plot_$date.plot

echo set output \"./$dirname/SpeedUp-OpenMP.png\" >> ./$dirname/plot_$date.plot
echo "plot x title 'Speedup max' with lines," \'./$dirname/res_$date.data\' "using 1:(\$2/\$3) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines," \'./$dirname/res_$date.data\' "using 1:(\$2/\$3/(\$1)) title 'Efficacité' with linespoints;" >> ./$dirname/plot_$date.plot

gnuplot ./$dirname/plot_$date.plot
