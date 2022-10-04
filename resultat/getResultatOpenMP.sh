#!/bin/bash

n_particule=1500
dt=1
for i in 1 2 3 4 5 6 7 8
do
  ../parallelOpenMP/nbody_brute_force $n_particule $dt $i 1
done
