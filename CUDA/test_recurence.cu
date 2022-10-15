#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#include "kernel.cu"

int *compteur;

void init(){
    compteur=0;
}

int main(int argc, char **argv){
    cudaMallocManaged(&compteur,sizeof(int));
    compteur=0;
    h_k_incremente<<<1,1>>>(compteur);
    cudaDeviceSynchronize();

    printf("%d\n",compteur);

    cudaFree(compteur);
    return 0;
}