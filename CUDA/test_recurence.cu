#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#include<math.h>

#define NSTEP 1000

int const graphEtage=2;
int const n = pow(4,graphEtage);
int n_node=0;

int *valeur;

__global__ void k_set(int *valeur,int index,int lap){
    valeur[index]=lap;
}


void init(){
    for (int i=0;i<graphEtage;i++){
        n_node+=pow(4,i);
    }
    printf("n_node %d\n",n_node);


}

void recursiveLaunch(int index,int lap){
    if(index>n_node-n){
        return;
    }
    k_set<<<1,1>>>(valeur,index,lap);
    recursiveLaunch(index+1,lap);
}

int main(int argc, char **argv){
    init();
    cudaMallocManaged(&valeur,n_node*sizeof(int));
    for (int i =0;i<n_node;i++){
        valeur[i]=-1;
    }

    bool graphCreated=false;
    cudaGraph_t graph;
    cudaGraphExec_t instance;
    cudaStream_t stream;

    for(int istep=0; istep<NSTEP; istep++){
        if(!graphCreated){
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

            recursiveLaunch(*valeur,istep);

            cudaStreamEndCapture(stream, &graph);
            cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);
            graphCreated=true;
        }
        cudaGraphLaunch(instance, stream);
        cudaStreamSynchronize(stream);
    }

    cudaFree(valeur);
    return 0;
}