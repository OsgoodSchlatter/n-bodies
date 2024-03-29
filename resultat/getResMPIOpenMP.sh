#!/bin/bash

#nbody_brute_force nbody_barnes_hut
algo=nbody_brute_force_async

MAX_PROCESS=6
THREAD_VARIABLE=0
MAX_THREAD=6

N_PARTICULE=3000
T_FINAL=1

date=$(date +"%d_%m_%y_%s")
hostsfile="./hosts"
dirname=MPI_OpenMP_$algo\_$date

seq_duration=0

#Build the structure of data
#./$date
#./$date/log
#./$date/res_$date.data

mkdir $dirname
mkdir $dirname/log

#mpirun de chauffe
OMP_NUM_THREADS=2 mpirun -n 2 -f $hostsfile ../MPI_OpenMP/$algo $N_PARTICULE $T_FINAL

#Param init
echo \#MAX_PROCESS MAX_THREAD N_PARTICULE T_FINAL THREAD_VARIABLE> ./$dirname/param.data
echo $MAX_PROCESS $MAX_THREAD $N_PARTICULE $T_FINAL $THREAD_VARIABLE> ./$dirname/param.data >&1
echo mpi brute force synchrone> ./$dirname/param.data >&1

echo \#n_process n_thread p t_seq t_parallel > ./$dirname/res_$date.data

#echo $n_process $MAX_THREAD
#../sequential/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_seq.log >&1
OMP_NUM_THREADS=1 mpirun -n 1 ../MPI_OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_seq.log >&1

seq_duration=$(cat ./$dirname/log/log_seq.log | grep "Simulation" | cut -d " " -f 3)

echo  1 1 1 $seq_duration $seq_duration >> ./$dirname/res_$date.data

for n_process in $(seq 1 $MAX_PROCESS)
do
  if [ $THREAD_VARIABLE == 0 ]
  then
    for n_thread in $(seq 1 $MAX_THREAD)
    do
      #echo $n_process $n_thread
      OMP_NUM_THREADS=$n_thread mpirun -n $n_process -f $hostsfile ../MPI_OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_$n_process\_$n_thread.log >&1

      duration=$(cat ./$dirname/log/log_$n_process\_$n_thread.log | grep "Simulation" | cut -d " " -f 3)

      echo $n_process $n_thread $(expr $n_process \* $n_thread) $seq_duration $duration >> ./$dirname/res_$date.data

      echo $n_process $n_thread

    done
  else
    #echo $n_process $MAX_THREAD
    OMP_NUM_THREADS=$MAX_THREAD mpirun -n $n_process -f $hostsfile ../MPI_OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_$n_process\_$MAX_THREAD.log >&1

    duration=$(cat ./$dirname/log/log_$n_process\_$MAX_THREAD.log | grep "Simulation" | cut -d " " -f 3)

    echo $n_process $MAX_THREAD $(expr $n_process \* $MAX_THREAD) $seq_duration $duration >> ./$dirname/res_$date.data
  fi
done

OMP_NUM_THREADS=$MAX_THREAD mpirun -n 10 -f $hostsfile ../MPI_OpenMP/$algo $N_PARTICULE $T_FINAL > ./$dirname/log/log_$n_process\_$MAX_THREAD.log >&1
duration=$(cat ./$dirname/log/log_$n_process\_$MAX_THREAD.log | grep "Simulation" | cut -d " " -f 3)
echo 10 $MAX_THREAD $(expr 10 \* $MAX_THREAD) $seq_duration $duration >> ./$dirname/res_$date.data

#AFFICHAGE PLOT
#$cp ./plot.plot ./$date/plot_$date.plot
echo set terminal png >> ./$dirname/plot_$date.plot
echo set ylabel \"speedup\" >> ./$dirname/plot_$date.plot
echo set xlabel \"number of cores\" >> ./$dirname/plot_$date.plot
echo set xtics 1 >> ./$dirname/plot_$date.plot
echo set title \"$algo Speedup with MPI + OpenMP execution\" >> ./$dirname/plot_$date.plot

echo set output \"./$dirname/SpeedUp-MPI_OpenMP.png\" >> ./$dirname/plot_$date.plot
echo "plot x title 'Speedup max' with lines," \'./$dirname/res_$date.data\' "using 3:(\$4/\$5) title 'SpeedUp' with linespoints, 1 title 'Efficacité max' with lines," \'./$dirname/res_$date.data\' "using 3:(\$4/\$5/(\$3)) title 'Efficacité' with linespoints;" >> ./$dirname/plot_$date.plot

gnuplot ./$dirname/plot_$date.plot
