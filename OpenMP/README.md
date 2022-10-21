#OpenMP

##Version disponible
nbody_brute_force

##Compilation
make  
Avec Eztrace : make CC="eztrace_cc gcc"

##Lancement
OMP_NUM_THREADS=[n_thread] ./nbody_brute_force [nparticules=1500] [TFINAL=1]

Avec Eztrace :
eztrace -t openmp ./nbody_brute_force [nparticules=1500] [TFINAL=1]

Lire trace otf2 :
vite nom_trace/eztrace_log.otf2
