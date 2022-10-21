# MPI

## Version disponible
nbody_brute_force nbody_barnes_hut nbody_brute_force_async  

## Compilation
make  
Avec Eztrace : make CC="eztrace_cc mpicc"  

## Lancement
mpirun -n [n_process] -f hosts ./nbody_brute_force [nparticules=1500] [TFINAL=1]  
mpirun -n [n_process] -f hosts ./nbody_barnes_hut [nparticules=1500] [TFINAL=1]  
mpirun -n [n_process] -f hosts ./nbody_brute_force_async [nparticules=1000] [TFINAL=1]  
mpirun -n [n_process] -f hosts ./nbody_barnes_hut_v2 [nparticules=1500] [TFINAL=1]  

__Avec Eztrace :__  
mpirun -n [n_process] -f hosts eztrace -t mpi ./nbody_brute_force [nparticules=1500] [TFINAL=1]  
mpirun -n [n_process] -f hosts eztrace -t mpi ./nbody_barnes_hut [nparticules=1500] [TFINAL=1]  
mpirun -n [n_process] -f hosts eztrace -t mpi ./nbody_brute_force_async [nparticules=1000] [TFINAL=1]  
mpirun -n [n_process] -f hosts eztrace -t mpi ./nbody_barnes_hut_v2 [nparticules=1500] [TFINAL=1]  

__Lire trace otf2 :__  
vite nom_trace/eztrace_log.otf2  
