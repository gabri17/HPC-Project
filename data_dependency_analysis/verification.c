#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv){

    if(argc != 3){
        return 1;
    }
    
    int processes = atoi(argv[1]);
    int k = atoi(argv[2]);

    if(processes < k){
        int temp = k;
        k = processes;
        processes = temp;
    }

    int indexOfProcess = processes;
    int left = k;
    for(int t = left; t > 0; t--){
        indexOfProcess = indexOfProcess - 1; 
        printf("indexOfProcesses %d at iteration %d\n", indexOfProcess, t);
    }

    printf("\n");

    for(int t = left; t > 0; t--){
        indexOfProcess = processes - left + t - 1; 
        printf("indexOfProcesses %d at iteration %d\n", indexOfProcess, t);
    }


    return 0;
}