#include <cuda.h>
#include <cuda_runtime.h>

__global__ void k_k_incremente(int* valeur){
    if(valeur[0]==5){
        return;
    }
    valeur+=1;
    k_k_incremente<<<1,1>>>(valeur);
}

__global__ void h_k_incremente(int* valeur){
    if(valeur[0]==5){
        return;
    }
    valeur+=1;
    k_k_incremente<<<1,1>>>(valeur);
}