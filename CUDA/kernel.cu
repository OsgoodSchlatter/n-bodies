#include <cuda.h>
#include <cuda_runtime.h>

__global__ void k_incremente(int* valeur){
    if(valeur[0]==5){
        return;
    }
    valeur+=1;
    k_incremente<<<1,1>>>(valeur);
}