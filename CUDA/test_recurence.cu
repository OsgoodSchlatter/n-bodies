#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

int *compteur;

__global__ void k_incremente(int* valeur){
    if(valeur==5){
        return 0;
    }
    valeur+=1;
    k_incremente(valeur);
}

void init(){
    compteur=0;
}

int main(int argc, char **argv){
    cudaMallocManaged(&compteur,sizeof(int));
    init()
    k_incremente<<<1,1>>>(compteur);
    cudaDeviceSynchronize();

    printf("%d\n",compteur);

    cudaFree(compteur);
    return 0;
}