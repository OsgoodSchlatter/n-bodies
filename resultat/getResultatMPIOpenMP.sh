#!/bin/bash

MAX_PROCESS=2
THREAD_VARIABLE=0
MAX_THREAD=2

N_PARTICULE=1500
T_FINAL=1

for n_process in $(seq 1 $MAX_PROCESS)
do
  if [ $THREAD_VARIABLE == 0 ]
  then
    for n_thread in $(seq 1 $MAX_THREAD)
    do
      echo $n_process $n_thread
      OMP_NUM_THREADS=$n_thread mpirun -n $n_process ../MPI_OpenMP/nbody_brute_force &N_PARTICULE &T_FINAL >> ./log/log_$n_process_$n_thread.data >&1
    done
  else
    echo $n_process $MAX_THREAD
  fi
done