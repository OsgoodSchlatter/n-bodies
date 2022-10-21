#Sequential

##Version disponible
nbody_brute_force nbody_barnes_hut  

##Compilation
make  
Avec Eztrace : make CC="eztrace_cc gcc"  

##Lancement
./nbody_brute_force [nparticules=1500] [TFINAL=1]  
./nbody_barnes_hut [nparticules=1000] [TFINAL=1]  

Avec Eztrace :  
eztrace ./nbody_brute_force [nparticules=1500] [TFINAL=1]  
eztrace ./nbody_barnes_hut [nparticules=1000] [TFINAL=1]  

Lire trace otf2 :  
vite nom_trace/eztrace_log.otf2  
