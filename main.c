//
// Created by Anastasios on 13/12/2018.
//
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "mpiFindMedian.h"


int main(int argc, char **argv) {

    int processID, noProcesses, size, loop, master; // Size = # of elems, partLength = length of partition in each process
    float distances[16] = {1987.625, 1226.345, 1990.625, 9256.975, 15846.75, 18799.0,
                           4545.875, 17029.75, 5136.125, 19161.75, 1099.545, 7943.25,
                           435035.0, 676696.5, 10271.04, 992468.5};
    size = 16;

    /*srand((processID + 1) * time(NULL));
    distances = (float *) malloc(size * sizeof(float));
    for (int i = 0; i < size; ++i) {
        distances[i] = ((float)rand())/1000 ;//- rand();
    }*/

    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);    /* get number of processes */

    loop = 0;
    noProcesses = noProcesses / 1; // /2
    size = size / 1; // /2
    master = processID / noProcesses; // Integer division <=> Floor
    MPI_Comm communicator;
    MPI_Comm_split(MPI_COMM_WORLD, master, 0, &communicator);
    mpiFindMedian(processID, master, noProcesses, size, distances, loop, communicator);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

}
