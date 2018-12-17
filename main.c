//
// Created by Anastasios on 13/12/2018.
//
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include "mpiFindMedian.h"

struct point {
    float x, y;
};

typedef struct point point;

int main(int argc, char **argv) {

    int processID, noProcesses;
    int totalSize, perProcessSize, loop, master, groups; // Size = # of elems

    float *distances;
    point *points;
    totalSize = 16;

    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);    /* get number of processes */
    perProcessSize = totalSize / noProcesses;

    // Each process read simultaneously its data, from the file
    FILE *fp;
    char buff[255], buff2[255];

    long int offset = ((long int) perProcessSize) * 2 * processID; // *2 because x and y are on different lines
    fp = fopen("cities.txt", "r");
    for (long int i = 0; i < offset; ++i) {
        fscanf(fp, "%s", buff); // Moving stream according to offset
    }

    points = (point *) malloc(perProcessSize * sizeof(point));

    for (int i = 0; i < totalSize / noProcesses; ++i) {
        fscanf(fp, "%s", buff);
        points[i].x = atof(buff);
        fscanf(fp, "%s", buff2);
        points[i].y = atof(buff2);
        printf("%ld:\t%f\t%f\n", offset / 2 + i + 1, points[i].x, points[i].y);
    }

    fclose(fp);

    //int loopend = log2(noProcesses);
    for (loop = 1; loop < 2; ++loop) {
        groups = (int) pow(2, loop);
        //printf("Groups = %d\n", groups);
        noProcesses = noProcesses / groups; // No of processes per group
        totalSize = totalSize / groups; // Size of the array for each group
        master = processID / noProcesses; // Integer division <=> Floor, for positives
        master = master * groups;// Calculate the master of the process
        //printf("My Master is = %d\n", master);

        MPI_Comm communicator[1];
        MPI_Comm_split(MPI_COMM_WORLD, master, 0, communicator);

        distances = (float *) malloc((totalSize / noProcesses) * sizeof(float));
        generateNumbers(distances, totalSize / noProcesses, processID);

        //float median = mpiFindMedian(processID, noProcesses, size, distances, communicator);
    }
    free(distances);

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Finalize();

}
