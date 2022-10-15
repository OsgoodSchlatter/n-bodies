#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#include "nbody.h"


int const *graphEtage;
int const *n;

node_t *root;


void init(){
    graphEtage=2;
    n = 4**(graphEtage)
    printf("%d\n",n);
}

int main(int argc, char **argv){
    //cudaMallocManaged(&compteur,sizeof(int));
    init();

    //cudaDeviceSynchronize();



    //cudaFree(compteur);
    return 0;
}