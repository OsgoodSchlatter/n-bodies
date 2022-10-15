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

__global__ void k_set(int *valeur,int index){
    valeur[index]=1;
}


void init(){
    for (int i=0;i<graphEtage;i++){
        n_node+=pow(4,i);
    }
    printf("n_node %d\n",n_node);


}

void recursiveLaunch(int actualGraphEtage){
    if(graphEtage==actualGraphEtage){
        return;
    }
    for (int i=0;i<4;i++){
        k_set<<<1,1>>>(valeur);
        recursiveLaunch(actualGraphEtage+1);
    }
}

int main(int argc, char **argv){
    init();
    cudaMallocManaged(&valeur,n_node*sizeof(int));
    for (int i =0;i<n_node;i++){
        valuer[0]=-1;
    }

    bool graphCreated=false;
    cudaGraph_t graph;
    cudaGraphExec_t instance;
    for(int istep=0; istep<NSTEP; istep++){
        if(!graphCreated){
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

            recusiveLaunch();

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