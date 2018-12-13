//
// Created by Anastasios on 13/12/2018.
//
#include <stdlib.h>
#include <mpi.h>
#include "mpiFindMedian.h"


int main(int argc, char **argv) {

    int processId, noProcesses, size; // Size = # of elems, partLength = length of partition in each process

    size = atoi(argv[1]); // Gets size from input
    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);    /* get number of processes */

    mpiFindMedian(processId,noProcesses,size);

    MPI_Finalize();

}
