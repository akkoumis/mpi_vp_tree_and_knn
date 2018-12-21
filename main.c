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

void map2dto1d(float **a, float *b, int dimension, int totalData) {
    for (int i = 0; i < totalData; ++i) {

    }
}

float calculateDistance(float *origin, float *end, int dimension) {
    float sum = 0;
    for (int i = 0; i < dimension; ++i) {
        sum += powf(origin[i] - end[i], 2);
    }
    return sqrt(sum);
}

void copyPoint(float *new, float *old, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        new[i] = old[i];
    }
}

int checkStatuses(int *statuses, int size) {
    for (int i = 0; i < size; ++i) {
        if (statuses[i] < 1)
            return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    MPI_Status mpiStat;

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

    offsets[0] = offsetof(struct point, x);
    offsets[1] = offsetof(struct point, y);

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
        //printf("pid = %d.\tCommunicator = %d.\n", pid, *communicator);


        //point vp;
        float vp[d];
        if (pid == 0) {
            for (int i = 0; i < d; ++i) {
                vp[i] = pointsCoordinates[0][i]; // TODO Select random instead of the first one
            }
        }
        MPI_Bcast(&vp, d, MPI_FLOAT, 0, *communicator); // printf("%d:\tx = %f\ty = %f\n", processID, vp.x, vp.y);

        for (int i = 0; i < perProcessSize; ++i) {
            distances[i] = calculateDistance(vp, pointsCoordinates[i], d);
            //printf("Distance to %d:\t%f\n", processID * perProcessSize + i + 1, distances[i]);
        }

        MPI_Barrier(*communicator);
        float median = mpiFindMedian(processID, noProcesses, perGroupSize, distances, communicator); // TODO use pid 1st
        float **lessEqual, **greater;
        lessEqual = (float **) malloc(perProcessSize * sizeof(float *));
        greater = (float **) malloc(perProcessSize * sizeof(float *));
        for (int k = 0; k < perProcessSize; ++k) {
            lessEqual[k] = (float *) malloc(d * sizeof(float));
            greater[k] = (float *) malloc(d * sizeof(float));
        }
        int countLessEqual = 0, countGreater = 0;
        int *countersLE, *countersG, *statuses;

        // Splits points and calculates the counters.
        for (int j = 0; j < perProcessSize; ++j) {
            if (distances[j] <= median) {
                copyPoint(lessEqual[countLessEqual], pointsCoordinates[j], d);
                countLessEqual++;
            } else {
                copyPoint(greater[countGreater], pointsCoordinates[j], d);
                countGreater++;
            }
        }

        if (pid == 0) {
            countersLE = (int *) malloc(noProcesses * sizeof(int));
            countersG = (int *) malloc(noProcesses * sizeof(int));
            statuses = (int *) malloc(noProcesses * sizeof(int));
            countersLE[0] = countLessEqual;
            countersG[0] = countGreater;
            statuses[0] = 0;
            for (int i = 1; i < noProcesses; ++i) {
                MPI_Recv(&countersLE[i], 1, MPI_INT, i, 3, *communicator, &mpiStat);
                MPI_Recv(&countersG[i], 1, MPI_INT, i, 4, *communicator, &mpiStat);
                statuses[i] = 0;
            }
        } else {
            MPI_Send(&countLessEqual, 1, MPI_INT, 0, 3, *communicator);
            MPI_Send(&countGreater, 1, MPI_INT, 0, 4, *communicator);
        }
        MPI_Barrier(*communicator);
        printf("ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", processID, countLessEqual, countGreater);
        MPI_Barrier(*communicator); // TODO Remove this. It's only for show.

        int status = 0, *directions;
        directions = (int *) malloc(3 * sizeof(int));

        if (pid == 0) {
            // MASTER: Swapping
            // Print state before swapping
            for (int i = 0; i < noProcesses; ++i) {
                printf("Master -> ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", i, countersLE[i], countersG[i]);
            }
            int indexLE = 0, indexG = noProcesses / 2;
            /*while (checkStatuses(statuses, noProcesses) == 0) { // While there are processes to be done.
                int *counterToBeZeroLE = (pid < noProcesses / 2) ? &countersG[indexLE] : &countersLE[indexLE];
                int *counterToBeMaxLE = (pid < noProcesses / 2) ? &countersLE[indexLE] : &countersG[indexG];
                int *counterToBeZeroG = (pid < noProcesses / 2) ? &countersG[indexG] : &countersLE[indexLE];
                int *counterToBeMaxG = (pid < noProcesses / 2) ? &countersLE[indexLE] : &countersG[indexG];
                if (*counterToBeZero == 0 && *counterToBeMax == perProcessSize){
                    ch
                }


                if (indexLE < noProcesses / 2) {
                    // Means that we are moving lessEqual elements to the first half
                    if (statuses[indexLE] == 2) {

                        indexLE++;
                    }
                } else {
                    // Means that we are moving greater elements to the second half
                    if (statuses[indexG] == 2)
                        indexG++;
                }
            }
            directions[0] = 3;
            directions[1] = 1;
            directions[2] = 1;
            countersG[3] -= 1;
            countersLE[1] += 1;
            MPI_Send(directions, 3, MPI_INT, 3, 6, *communicator);
            MPI_Send(directions, 3, MPI_INT, 1, 6, *communicator);*/

        } else {
            // SLAVE: Swapping
            int *counterToBeZero = (pid < noProcesses / 2) ? &countGreater : &countLessEqual;
            int *counterToBeMax = (pid < noProcesses / 2) ? &countLessEqual : &countGreater;
            while (1) {
                //MPI_Recv(&status, 1, MPI_INT, 0, 5, *communicator, &mpiStat);
                if (*counterToBeZero == 0 && *counterToBeMax == perProcessSize)
                    break;
                MPI_Recv(directions, 3, MPI_INT, 0, 6, *communicator, &mpiStat);
                int totalData = directions[2] * d;
                if (directions[0] == pid) {
                    // If this SLAVE must SEND
                    float **dataToSend = (pid < noProcesses / 2) ? greater : lessEqual;
                    int *counterToModify = (pid < noProcesses / 2) ? &countGreater : &countLessEqual;
                    *counterToModify -= directions[2];
                    float *arrayToSend = (float *) malloc(totalData * sizeof(float));
                    for (int i = 0; i < directions[2]; ++i) {
                        for (int j = 0; j < d; ++j) {
                            arrayToSend[i * d + j] = dataToSend[*counterToModify + i][j];
                        }
                    }
                    MPI_Send(arrayToSend, totalData, MPI_FLOAT, directions[1], 7, *communicator);
                } else {
                    // If this SLAVE must RECEIVE
                    float **dataToReceive = (pid < noProcesses / 2) ? lessEqual : greater;
                    int *counterToModify = (pid < noProcesses / 2) ? &countLessEqual : &countGreater;
                    float *arrayToReceive = (float *) malloc(totalData * sizeof(float));
                    MPI_Recv(arrayToReceive, totalData, MPI_FLOAT, directions[0], 7, *communicator, &mpiStat);
                    for (int i = 0; i < directions[2]; ++i) {
                        for (int j = 0; j < d; ++j) {
                            dataToReceive[*counterToModify + i][j] = arrayToReceive[i * d + j];
                        }
                    }
                    *counterToModify += directions[2];
                }
                printf("ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", processID, countLessEqual, countGreater);
            }
        }


    }

    // TODO More frees of pointers with malloc
    free(distances);
    for (int k = 0; k < perProcessSize; ++k) {
        free(pointsCoordinates[k]);
    }
    free(pointsCoordinates);
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Type_free(&mpi_point_type);
    MPI_Finalize();

}
