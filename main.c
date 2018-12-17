//
// Created by Anastasios on 13/12/2018.
//
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
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

    /* create a type for struct point */
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_point_type;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(struct point, x);
    offsets[1] = offsetof(struct point, y);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);


    // Each process read simultaneously its data, from the file
    FILE *fp;
    char buff[255], buff2[255];

    long int offset = ((long int) perProcessSize) * 2 * processID; // *2 because x and y are on different lines
    fp = fopen("cities.txt", "r");
    for (long int i = 0; i < offset; ++i) {
        fscanf(fp, "%s", buff); // Moving stream according to offset
    }

    points = (point *) malloc(perProcessSize * sizeof(point));
    distances = (float *) malloc(perProcessSize * sizeof(float));

    for (int i = 0; i < totalSize / noProcesses; ++i) {
        fscanf(fp, "%s", buff);
        points[i].x = atof(buff);
        fscanf(fp, "%s", buff2);
        points[i].y = atof(buff2);
        printf("%ld:\t%f\t%f\n", offset / 2 + i + 1, points[i].x, points[i].y);
    }
    fclose(fp);


    //int loopend = log2(noProcesses);
    for (loop = 0; loop < 1; ++loop) {
        groups = (int) pow(2, loop);
        //printf("Groups = %d\n", groups);
        noProcesses = noProcesses / groups; // No of processes per group
        totalSize = totalSize / groups; // Size of the array for each group
        master = processID / noProcesses; // Integer division <=> Floor, for positives
        master = master * groups;// Calculate the master of the process
        //printf("My Master is = %d\n", master);

        MPI_Comm communicator[1];
        MPI_Comm_split(MPI_COMM_WORLD, master, 0, communicator);

        point vp;
        if (processID == 0) {
            vp.x = points[0].x;
            vp.y = points[0].y;
        }
        MPI_Bcast(&vp, 1, mpi_point_type, 0, *communicator);
        printf("%d:\tx = %f\ty = %f\n", processID, vp.x, vp.y);

        distances = (float *) malloc((totalSize / noProcesses) * sizeof(float));
        generateNumbers(distances, totalSize / noProcesses, processID);

        //float median = mpiFindMedian(processID, noProcesses, size, distances, communicator);
    }
    free(distances);

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Type_free(&mpi_point_type);
    MPI_Finalize();

}
