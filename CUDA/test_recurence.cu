#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#include "nbody.h"


int const graphEtage=2;
int const n = 4**(graphEtage);

node_t *root;


void init(){
    printf("%d\n",n);
}

int main(int argc, char **argv){
    //cudaMallocManaged(&compteur,sizeof(int));
    init();

    //cudaDeviceSynchronize();



    //cudaFree(compteur);
    return 0;
}