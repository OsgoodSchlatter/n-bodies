# MPI_OpenMP

## Version disponible
nbody_brute_force nbody_barnes_hut nbody_brute_force_async  

## Compilation
make  
_Avec Eztrace_ : make CC="eztrace_cc mpicc"  

## Lancement
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts ./nbody_brute_force [nparticules=1500] [TFINAL=1]  
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts ./nbody_barnes_hut [nparticules=1500] [TFINAL=1]  
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts ./nbody_brute_force_async [nparticules=1000] [TFINAL=1]  

_Avec Eztrace_ :  
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts eztrace -t "mpi openmp" ./nbody_brute_force [nparticules=1500] [TFINAL=1]  
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts eztrace -t "mpi openmp" ./nbody_barnes_hut [nparticules=1500] [TFINAL=1]  
OMP_NUM_THREADS=[n_thread] mpirun -n [n_process] -f hosts eztrace -t "mpi openmp" ./nbody_brute_force_async [nparticules=1000] [TFINAL=1]  

_Lire trace otf2_ :  
vite nom_trace/eztrace_log.otf2  
