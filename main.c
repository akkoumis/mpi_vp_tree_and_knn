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

int main(int argc, char **argv) {

    int processID, noProcesses;
    int size, loop, master, groups; // Size = # of elems, partLength = length of partition in each process

    float *distances;
    struct point *points;
    size = 16;

    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);    /* get number of processes */

    // Each process read simultaneously its data, from the file
    FILE *fp;
    char buff[255], buff2[255];
    float coordinate, coordinate2;
    long int offset = (((long int) size) / noProcesses) * 2 * processID; // *2 because x and y are on different lines

    fp = fopen("cities.txt", "r");
    //fseek(fp, offset, SEEK_CUR);

    for (long int i = 0; i < offset; ++i) {
        fscanf(fp, "%s", buff);
    }

    for (int i = 0; i < size / noProcesses; ++i) {
        fscanf(fp, "%s", buff);
        coordinate = atof(buff);
        fscanf(fp, "%s", buff2);
        coordinate2 = atof(buff2);
        printf("%d:\t%f\t%f\n", offset / 2 + i + 1, coordinate, coordinate2);
    }

    fclose(fp);

    //int loopend = log2(noProcesses);
    for (loop = 1; loop < 2; ++loop) {
        groups = (int) pow(2, loop);
        //printf("Groups = %d\n", groups);
        noProcesses = noProcesses / groups; // No of processes per group
        size = size / groups; // Size of the array for each group
        master = processID / noProcesses; // Integer division <=> Floor, for positives
        master = master * groups;// Calculate the master of the process
        //printf("My Master is = %d\n", master);

        MPI_Comm communicator[1];
        MPI_Comm_split(MPI_COMM_WORLD, master, 0, communicator);

        distances = (float *) malloc((size / noProcesses) * sizeof(float));
        generateNumbers(distances, size / noProcesses, processID);

        //float median = mpiFindMedian(processID, noProcesses, size, distances, communicator);
    }
    free(distances);

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Finalize();

}
