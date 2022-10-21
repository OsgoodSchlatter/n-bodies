# OpenMP

## Version disponible
nbody_brute_force

## Compilation
make  
_Avec Eztrace_ : make CC="eztrace_cc gcc"

## Lancement
OMP_NUM_THREADS=[n_thread] ./nbody_brute_force [nparticules=1500] [TFINAL=1]

_Avec Eztrace_ :
eztrace -t openmp ./nbody_brute_force [nparticules=1500] [TFINAL=1]

_Lire trace otf2_ :
vite nom_trace/eztrace_log.otf2
