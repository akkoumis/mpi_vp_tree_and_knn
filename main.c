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

float calculateDistance(float *origin, float *end, int dimension) {
    float sum = 0;
    for (int i = 0; i < dimension; ++i) {
        sum += powf(origin[i] - end[i], 2);
    }
    return sqrt(sum);
}

int main(int argc, char **argv) {
    MPI_Status Stat;

    int d = 2;
    int processID, pid, noProcesses, noTotalProcesses;
    int totalSize, perProcessSize, loop, master, groups, perGroupSize; // Size = # of elems

    float *distances, **pointsCoordinates;
    point *points;
    totalSize = 16;

    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noTotalProcesses);    /* get total number of processes */
    perProcessSize = totalSize / noTotalProcesses;

    /* create a type for struct point */
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_point_type;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(
    struct point, x);
    offsets[1] = offsetof(
    struct point, y);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);


    // Each process read simultaneously its data, from the file
    FILE *fp;
    char buff[255], buff2[255];

    long int offset = ((long int) perProcessSize) * 2 * processID; // *2 -> because x and y are on different lines
    fp = fopen("cities.txt", "r");
    for (long int i = 0; i < offset; ++i) {
        fscanf(fp, "%s", buff); // Moving stream according to offset
    }

    //points = (point *) malloc(perProcessSize * sizeof(point));
    distances = (float *) malloc(perProcessSize * sizeof(float));
    pointsCoordinates = (float **) malloc(perProcessSize * sizeof(float *));
    for (int k = 0; k < perProcessSize; ++k) {
        pointsCoordinates[k] = (float *) malloc(d * sizeof(float));
    }

    for (int i = 0; i < perProcessSize; ++i) {
        for (int j = 0; j < d; ++j) {
            fscanf(fp, "%s", buff2);
            pointsCoordinates[i][j] = atof(buff2);
        }
        //printf("%ld:\t%f\t%f\n", offset / 2 + i + 1, points[i].x, points[i].y);
    }
    fclose(fp);

    MPI_Comm communicator[1];

    //int loopend = log2(noProcesses);
    for (loop = 0; loop <= 0; ++loop) {
        groups = (int) pow(2, loop); //number of groups
        //printf("Groups = %d\n", groups);
        noProcesses = noTotalProcesses / groups; // No of processes per group
        perGroupSize = totalSize / groups; // Size of the array for each group
        master = processID / noProcesses; // Integer division <=> Floor, for positives
        master = master * noProcesses;// Calculate the master of the process
        //printf("My Master is = %d\n", master);

        //printf("Pid = %d.\tCommunicator = %d.\n", processID, *communicator);
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_split(MPI_COMM_WORLD, master, 0, communicator);
        MPI_Comm_rank(*communicator, &pid);
        //if(processID==1)
        //printf("Pid = %d.\tCommunicator = %d.\n", processID, *communicator);


        //point vp;
        float vp[d];
        if (pid == 0) {
            for (int i = 0; i < d; ++i) {
                vp[i] = pointsCoordinates[0][i]; // TODO Select random instead of the first one
            }
        }
        MPI_Bcast(&vp, d, MPI_FLOAT, 0, *communicator); // printf("%d:\tx = %f\ty = %f\n", processID, vp.x, vp.y);

        distances = (float *) malloc((perProcessSize) * sizeof(float));
        for (int i = 0; i < perProcessSize; ++i) {
            distances[i] = calculateDistance(vp, pointsCoordinates[i], d);
            //printf("Distance to %d:\t%f\n", processID * perProcessSize + i + 1, distances[i]);
        }

        MPI_Barrier(*communicator);
        float median = mpiFindMedian(processID, noProcesses, perGroupSize, distances, communicator); // TODO use pid 1st

        float *lessEqual, *greater;
        lessEqual = (float *) malloc(perGroupSize * sizeof(float));
        greater = (float *) malloc(perGroupSize * sizeof(float));
        int countLessEqual = 0, countGreater = 0;
        int *countersLE, *countersG;

        // Splits points and calculates the counters.
        for (int j = 0; j < perProcessSize; ++j) {
            if (distances[j] <= median) {
                lessEqual[countLessEqual] = distances[j];
                countLessEqual++;
            } else {
                greater[countGreater] = distances[j];
                countGreater++;
            }
        }

        if (pid == 0) {
            countersLE = (int *) malloc(noProcesses * sizeof(int));
            countersG = (int *) malloc(noProcesses * sizeof(int));
            countersLE[0] = countLessEqual;
            countersG[0] = countGreater;
            for (int i = 1; i < noProcesses; ++i) {
                MPI_Recv(&countersLE[i], 1, MPI_INT, i, 3, *communicator, &Stat);
                MPI_Recv(&countersG[i], 1, MPI_INT, i, 4, *communicator, &Stat);
            }
        } else {
            MPI_Send(&countLessEqual, 1, MPI_INT, 0, 3, *communicator);
            MPI_Send(&countGreater, 1, MPI_INT, 0, 4, *communicator);
        }
        MPI_Barrier(*communicator);
        printf("ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", processID, countLessEqual, countGreater);

        if (pid == 0) {
            for (int i = 0; i < noProcesses; ++i) {
                printf("Master -> ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", i, countersLE[i], countersG[i]);
            }
        }


    }
    free(distances);

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Type_free(&mpi_point_type);
    MPI_Finalize();

}
