#!/bin/bash

MAX_PROCESS=4
THREAD_VARIABLE=1
MAX_THREAD=2

N_PARTICULE=1500
T_FINAL=1
date=$(date +"%d_%m_%y_%s")

seq_duration=0

#Build the structure of data
#./$date
#./$date/log
#./$date/res_$date.data

mkdir $date
mkdir $date/log

echo \#n_process n_thread t_seq t_parallel > ./$date/res_$date.data

for n_process in $(seq 1 $MAX_PROCESS)
do
  if [ $THREAD_VARIABLE == 0 ]
  then
    for n_thread in $(seq 1 $MAX_THREAD)
    do
      #echo $n_process $n_thread
      OMP_NUM_THREADS=$n_thread mpirun -n $n_process ../MPI_OpenMP/nbody_brute_force $N_PARTICULE $T_FINAL > ./$date/log/log_$n_process\_$n_thread.log >&1

      duration=$(cat ./log/log_$n_process\_$n_thread.log | grep "Simulation" | cut -d " " -f 3)

      if [ $seq_duration == 0 ]
      then
        seq_duration=$duration
      fi

      echo $n_process $n_thread $seq_duration $duration >> ./res_$date.data

    done
  else
    #echo $n_process $MAX_THREAD
    OMP_NUM_THREADS=$MAX_THREAD mpirun -n $n_process ../MPI_OpenMP/nbody_brute_force $N_PARTICULE $T_FINAL > ./$date/log/log_$n_process\_$MAX_THREAD.log >&1

    duration=$(cat ./$date/log/log_$n_process\_$MAX_THREAD.log | grep "Simulation" | cut -d " " -f 3)

    if [ $seq_duration == 0 ]
    then
      seq_duration=$duration
    fi

    echo $n_process $MAX_THREAD $seq_duration $duration >> ./$date/res_$date.data
  fi
done

#AFFICHAGE PLOT
#$cp ./plot.plot ./$date/plot_$date.plot
echo set terminal png >> ./$date/plot_$date.plot
echo set ylabel \"speedup\" >> ./$date/plot_$date.plot
echo set xlabel \"number of cores\" >> ./$date/plot_$date.plot
echo set xtics 1 >> ./$date/plot_$date.plot
echo set title \"NBody_Brute_Force Speedup for MPI + OpenMP execution\" >> ./$date/plot_$date.plot

echo set output \"./$date/SpeedUp-MPI_OpenMP.png\" >> ./$date/plot_$date.plot
echo "plot x title 'Speedup max' with lines," \'./$date/res_$date.data\' "using 1:(\$3/\$4) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines," \'./$date/res_$date.data\' "using 1:(\$3/\$4/(\$1)) title 'Efficacité' with linespoints;" >> ./$date/plot_$date.plot

gnuplot ./$date/plot_$date.plot
