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

struct node {
    float *vp;
    float radius;
};

typedef struct node node;


int findMin(int a, int b) {
    if (a < b)
        return a;
    return b;
}

float calculateDistance(const float *origin, const float *end, int dimension) {
    float sum = 0;
    for (int i = 0; i < dimension; ++i) {
        sum += powf(origin[i] - end[i], 2);
    }
    float res = sqrt(sum);
    //printf("Distance to (%f,%f):\t%f\n", end[0], end[1], res);
    return res;
}

void copyPoint(float *new, float *old, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        new[i] = old[i];
    }
}

int checkStatuses(const int *statuses, int size) {
    for (int i = 0; i < size; ++i) {
        if (statuses[i] < 1)
            return 0;
    }
    return 1;
}

void finishTree(int processID, node *tree, int index, float **points, int noOfPoints, int dimension, MPI_Comm *communicator) {
    if (noOfPoints == 1) {
        tree[index].radius = -1;                // Insert node's radius
        tree[index].vp = (float *) malloc(dimension * sizeof(float));
        copyPoint(tree[index].vp, points[0], dimension); // Insert the node's vp
        //printf("tree[%d].vp = (%f,%f)\n", index, tree[index].vp[0], tree[index].vp[1]);
        return;
    } else {
        float *distances, **lessEqual, **greater;
        int countLessEqual = 0, countGreater = 0;
        distances = (float *) malloc(noOfPoints * sizeof(float)); // Allocate memory for distances
        lessEqual = (float **) malloc((noOfPoints / 2) * sizeof(float *));
        greater = (float **) malloc((noOfPoints / 2) * sizeof(float *));
        for (int k = 0; k < noOfPoints / 2; ++k) {
            lessEqual[k] = (float *) malloc(dimension * sizeof(float));
            greater[k] = (float *) malloc(dimension * sizeof(float));
        }

        float *vp = points[0]; // Select VP TODO: Select randomly
        for (int i = 0; i < noOfPoints; ++i) {
            distances[i] = calculateDistance(vp, points[i], dimension); // Calculate distances
            //printf("processID %d\tdistances[%d]=%f\n", processID, i, distances[i]);
        }

        float median = mpiFindMedian(0, 1, noOfPoints, distances, communicator);
        tree[index].vp = (float *) malloc(dimension * sizeof(float));
        tree[index].radius = median; // Insert node's radius
        copyPoint(tree[index].vp, vp, dimension); // Insert the node's vp
        for (int j = 0; j < noOfPoints; ++j) { // Split points according to distance from median
            if (distances[j] <= median) {
                copyPoint(lessEqual[countLessEqual], points[j], dimension);
                countLessEqual++;
            } else {
                copyPoint(greater[countGreater], points[j], dimension);
                countGreater++;
            }
        }

        if (countGreater != (noOfPoints / 2) || countLessEqual != (noOfPoints / 2))
            printf("ERROR splitting finishTree() in pid = %d cle=%d cg=%d\n", processID, countLessEqual, countGreater);

        finishTree(processID, tree, 2 * index + 1, lessEqual, noOfPoints / 2, dimension, communicator);
        finishTree(processID, tree, 2 * index + 2, greater, noOfPoints / 2, dimension, communicator);

        free(lessEqual);
        free(greater);
        free(distances);
    }
}

void shareTree(int myProcessID, int sendersProcessID, node *tree, int index, int treeSize, int dimension) {
    MPI_Bcast(&(tree[index].radius), 1, MPI_FLOAT, sendersProcessID, MPI_COMM_WORLD);
    if (myProcessID != sendersProcessID)
        tree[index].vp = (float *) malloc(dimension * sizeof(float));
    MPI_Bcast(tree[index].vp, dimension, MPI_FLOAT, sendersProcessID, MPI_COMM_WORLD);
    if ((index * 2 + 2) < treeSize) {
        shareTree(myProcessID, sendersProcessID, tree, index * 2 + 1, treeSize, dimension);
        shareTree(myProcessID, sendersProcessID, tree, index * 2 + 2, treeSize, dimension);
    }
}

int main(int argc, char **argv) {
    MPI_Status mpiStat;

    int d = 2;
    int processID, pid, noProcesses, noTotalProcesses;
    int totalSize, perProcessSize, loop, master, groups, perGroupSize; // Size = # of elems

    float *distances, **pointsCoordinates, **lessEqual, **greater;;
    totalSize = 16;

    MPI_Init(&argc, &argv);    /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);    /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &noTotalProcesses);    /* get total number of processes */
    perProcessSize = totalSize / noTotalProcesses;

    /* create a type for struct node */
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_LONG_INT, MPI_FLOAT};
    MPI_Datatype mpi_node_type;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(struct node, vp);
    offsets[1] = offsetof(struct node, radius);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_node_type);
    MPI_Type_commit(&mpi_node_type);

    //printf("%d\t%d\n", sizeof(node), sizeof(long int));

    node *tree = (node *) malloc((2 * totalSize - 1) * sizeof(node));
    //printf("tree size = %d\n", (2 * totalSize - 1));

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
    float *tempDistances = (float *) malloc(perProcessSize * sizeof(float));
    pointsCoordinates = (float **) malloc(perProcessSize * sizeof(float *));
    lessEqual = (float **) malloc(perProcessSize * sizeof(float *));
    greater = (float **) malloc(perProcessSize * sizeof(float *));
    for (int k = 0; k < perProcessSize; ++k) {
        pointsCoordinates[k] = (float *) malloc(d * sizeof(float));
        lessEqual[k] = (float *) malloc(d * sizeof(float));
        greater[k] = (float *) malloc(d * sizeof(float));
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

    int loopEnd = (int) log2(noTotalProcesses);
    //printf("loopEnd = %d\n",loopEnd);
    for (loop = 0; loop < loopEnd; ++loop) {
        groups = (int) pow(2, loop); // Number of groups of processes currently working
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
            tempDistances[i] = distances[i];
            //printf("Distance to %d-%d (%f,%f):\t%f\n", processID, i, pointsCoordinates[i][0], pointsCoordinates[i][1],
            //       distances[i]);
        }


        MPI_Barrier(*communicator);
        float median = mpiFindMedian(processID, noProcesses, perGroupSize, tempDistances, communicator); // TODO use pid 1st

        int countLessEqual = 0, countGreater = 0;
        int *countersLE, *countersG, *statuses;

        //MPI_Barrier(*communicator);

        // Splits points and calculates the counters.
        for (int j = 0; j < perProcessSize; j++) {
            //printf("Process %d-%d: %f\t%f\tdist=%f\n", pid, j, pointsCoordinates[j][0], pointsCoordinates[j][1], distances[j]);
            if (distances[j] <= median) {
                copyPoint(lessEqual[countLessEqual], pointsCoordinates[j], d);
                //printf("Process %d-%d: lessEqual[%d] = %f\t%f\tdist=%f\n", pid, j, countLessEqual, lessEqual[countLessEqual][0],
                //       lessEqual[countLessEqual][1], distances[j]);
                countLessEqual++;
            } else {
                copyPoint(greater[countGreater], pointsCoordinates[j], d);
                //printf("Process %d-%d: greater[%d] = %f\t%f\tdist=%f\n", pid, j, countGreater, greater[countGreater][0],
                //       greater[countGreater][1], distances[j]);
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
        //printf("(BEFORE) ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", processID, countLessEqual, countGreater);
        //MPI_Barrier(*communicator); // TODO Remove this. It's only for show.

        int status = 0, directions[3];
        //directions = (int *) malloc(3 * sizeof(int));



        if (pid == 0) {
            // MASTER: Swapping
            // Print state before swapping
            for (int i = 0; i < noProcesses; ++i) {
                //printf("(BEFORE) Master -> ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", i, countersLE[i], countersG[i]);
            }
            int index1LE = 0, index1G = 0, index2LE = noProcesses / 2, index2G = noProcesses / 2;
            while (checkStatuses(statuses, noProcesses) == 0) { // While there are processes to be done.
                // Index1
                if (countersLE[index1LE] == perProcessSize) {
                    if (countersG[index1LE] == 0) {
                        statuses[index1LE] = 1;
                        //printf("pid(LE) = %d done!\n",index1LE);
                        MPI_Send(&statuses[index1LE], 1, MPI_INT, index1LE, 5, *communicator);
                    }
                    index1LE++;
                }
                if (countersG[index1G] == 0) {
                    if (countersLE[index1G] == perProcessSize) {
                        statuses[index1G] = 1;
                        //printf("pid(G) = %d done!\n",index1G);
                        MPI_Send(&statuses[index1G], 1, MPI_INT, index1G, 5, *communicator);
                    }
                    index1G++;
                }
                // Index2
                if (countersLE[index2LE] == 0) {
                    if (countersG[index2LE] == perProcessSize) {
                        statuses[index2LE] = 1;
                        //printf("pid(LE) = %d done!\n",index2LE);
                        MPI_Send(&statuses[index2LE], 1, MPI_INT, index2LE, 5, *communicator);
                    }
                    index2LE++;
                }
                if (countersG[index2G] == perProcessSize) {
                    if (countersLE[index2G] == 0) {
                        statuses[index2G] = 1;
                        //printf("pid(G) = %d done!\n",index2G);
                        MPI_Send(&statuses[index2G], 1, MPI_INT, index2G, 5, *communicator);
                    }
                    index2G++;
                }
                // Now we have indices that point to processes with data to be moved.

                if (index1LE < noProcesses / 2) { // We have to send LE to 1st subgroup
                    // At first SEND from index2LE and RECEIVE to index1LE
                    int LEToFull = perProcessSize - countersLE[index1LE];
                    int GElems = countersG[index2LE];
                    int noOfElements = findMin(LEToFull, GElems);
                    directions[0] = index2LE;
                    directions[1] = index1LE;
                    directions[2] = noOfElements;
                    status = 0;
                    MPI_Send(&status, 1, MPI_INT, directions[0], 5, *communicator);
                    MPI_Send(directions, 3, MPI_INT, directions[0], 6, *communicator); // Send directions to sender.
                    if (index1LE == 0) {
                        // MASTER must RECEIVE
                        // We have to receive as the slaves do
                        int totalData = directions[2] * d;
                        float **dataToModify = lessEqual;
                        int *counterToModify = &countLessEqual;
                        float *arrayToReceive = (float *) malloc(totalData * sizeof(float));
                        MPI_Recv(arrayToReceive, totalData, MPI_FLOAT, directions[0], 7, *communicator, &mpiStat);
                        for (int i = 0; i < directions[2]; ++i) {
                            for (int j = 0; j < d; ++j) {
                                dataToModify[*counterToModify + i][j] = arrayToReceive[i * d + j];
                                //printf("DATA RECV: pid=%d\t[%d][%d]\t%f\n", pid, *counterToModify + i, j, arrayToReceive[i * d + j]);
                            }
                        }
                        free(arrayToReceive);
                        *counterToModify += directions[2];
                    } else {
                        // We just send directions to the RECEIVING slave.
                        MPI_Send(&status, 1, MPI_INT, directions[1], 5, *communicator);
                        MPI_Send(directions, 3, MPI_INT, directions[1], 6, *communicator);
                    }
                    countersLE[index2LE] -= directions[2];
                    countersLE[index1LE] += directions[2];
                }

                if (index1G < noProcesses / 2) { // It means that we are not done
                    // Secondly SEND from index1G and RECEIVE to index2G
                    int GToFull = perProcessSize - countersG[index2G];
                    int LEElems = countersG[index1G];
                    int noOfElements = findMin(GToFull, LEElems);
                    directions[0] = index1G;
                    directions[1] = index2G;
                    directions[2] = noOfElements;
                    status = 0;
                    MPI_Send(&status, 1, MPI_INT, directions[1], 5, *communicator);
                    MPI_Send(directions, 3, MPI_INT, directions[1], 6, *communicator); // Send directions to receiver.
                    if (index1G == 0) {
                        // MASTER must SEND
                        int totalData = directions[2] * d;
                        float **dataToModify = greater;
                        int *counterToModify = &countGreater;
                        *counterToModify -= directions[2];
                        float *arrayToSend = (float *) malloc(totalData * sizeof(float));
                        for (int i = 0; i < directions[2]; ++i) {
                            for (int j = 0; j < d; ++j) {
                                arrayToSend[i * d + j] = dataToModify[*counterToModify + i][j];
                            }
                        }
                        MPI_Send(arrayToSend, totalData, MPI_FLOAT, directions[1], 7, *communicator);
                        for (int k = 0; k < directions[2]; ++k) {
                            //printf("SENT: pid=%d\t%f\t%f\n", pid, arrayToSend[k * d + 0], arrayToSend[k * d + 1]); // CHECKED
                        }
                        free(arrayToSend);
                    } else {
                        // We just send directions to the SENDING slave.
                        MPI_Send(&status, 1, MPI_INT, directions[0], 5, *communicator);
                        MPI_Send(directions, 3, MPI_INT, directions[0], 6, *communicator);
                    }
                    countersG[index1G] -= directions[2];
                    countersG[index2G] += directions[2];
                }
            }
            for (int i = 0; i < noProcesses; ++i) {
                //printf("(AFTER) Master -> ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", i, countersLE[i],countersG[i]);

            }
            /*directions[0] = 3;
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
                MPI_Recv(&status, 1, MPI_INT, 0, 5, *communicator, &mpiStat);
                if (status == 1)
                    break;
                if (*counterToBeZero == 0 && *counterToBeMax == perProcessSize)
                    break;
                MPI_Recv(directions, 3, MPI_INT, 0, 6, *communicator, &mpiStat); // Receive directions
                int totalData = directions[2] * d;

                if (directions[0] == pid) {
                    // If this SLAVE must SEND
                    float **dataToModify = (pid < noProcesses / 2) ? greater : lessEqual;
                    int *counterToModify = (pid < noProcesses / 2) ? &countGreater : &countLessEqual;
                    *counterToModify -= directions[2];
                    float *arrayToSend = (float *) malloc(totalData * sizeof(float));
                    //printf("SENDING:   pid=%d\t%d->%d\t%d elements from index %d\n", pid, directions[0], directions[1], directions[2],
                    //       *counterToModify);
                    for (int i = 0; i < directions[2]; ++i) {
                        for (int j = 0; j < d; ++j) {
                            arrayToSend[i * d + j] = dataToModify[*counterToModify + i][j];
                        }
                    }
                    MPI_Send(arrayToSend, totalData, MPI_FLOAT, directions[1], 7, *communicator);
                    for (int k = 0; k < directions[2]; ++k) {
                        //printf("SENT: pid=%d\t%f\t%f\n", pid, arrayToSend[k * d + 0], arrayToSend[k * d + 1]); // CHECKED
                    }
                    free(arrayToSend);
                } else {
                    // If this SLAVE must RECEIVE
                    float **dataToModify = (pid < noProcesses / 2) ? lessEqual : greater;
                    int *counterToModify = (pid < noProcesses / 2) ? &countLessEqual : &countGreater;
                    float *arrayToReceive = (float *) malloc(totalData * sizeof(float));
                    //printf("RECEIVING: pid=%d\t%d->%d\t%d elements to index %d\n", pid, directions[0], directions[1], directions[2],
                    //       *counterToModify);
                    MPI_Recv(arrayToReceive, totalData, MPI_FLOAT, directions[0], 7, *communicator, &mpiStat);
                    for (int k = 0; k < directions[2]; ++k) {
                        //printf("RECV: pid=%d\t%f\t%f\n",pid,arrayToReceive[k*d+0],arrayToReceive[k*d+1]);
                    }
                    for (int i = 0; i < directions[2]; ++i) {
                        for (int j = 0; j < d; ++j) {
                            dataToModify[*counterToModify + i][j] = arrayToReceive[i * d + j];
                            //printf("DATA RECV: pid=%d\t[%d][%d]\t%f\n", pid, *counterToModify + i, j, dataToModify[*counterToModify +
                            //i][j]);
                        }
                    }
                    free(arrayToReceive);
                    *counterToModify += directions[2];
                    if (*counterToModify == perProcessSize) {
                        for (int i = 0; i < perProcessSize; ++i) {
                            //printf("DATA END: pid=%d\t[%d]\t%f\t%f\n", pid, i, dataToModify[i][0], dataToModify[i][1]);
                        }
                    }
                }
            }
            if (*counterToBeZero > 0) {// It means we have some extra medians that shall remain here (2nd subgroup)
                //TODO
                printf("ERROR: Multiple medians\n");
                for (int i = 0; i < *counterToBeZero; ++i) {
                    for (int j = 0; j < d; ++j) {
                        greater[*counterToBeMax + i][j] = lessEqual[i][j]; // Move points from lessEqual to greater
                    }
                }
                *counterToBeMax += *counterToBeZero;
                *counterToBeZero = 0;
            }
        }
        float **dataToFetch = (pid < noProcesses / 2) ? lessEqual : greater;
        for (int l = 0; l < perProcessSize; ++l) {
            for (int i = 0; i < d; ++i) {
                pointsCoordinates[l][i] = dataToFetch[l][i]; // Save the data for the next iteration of the loop
            }
        }
        // Calculates the counters for verification.
        int tempCountLessEqual = 0, tempCountGreater = 0;
        for (int j = 0; j < perProcessSize; ++j) {
            if ((calculateDistance(vp, lessEqual[j], d) <= median) && pid < noProcesses / 2 && j < countLessEqual) {
                //printf("(AFTER) Process %d: lessEqual[%d] = %f\t%f\n", pid, tempCountLessEqual, pointsCoordinates[tempCountLessEqual][0],
                //       pointsCoordinates[tempCountLessEqual][1]);
                tempCountLessEqual++;
            }
            if ((calculateDistance(vp, greater[j], d) >= median) && pid >= noProcesses / 2 && j < countGreater) {
                //printf("(AFTER) Process %d: greater[%d] = %f\t%f\n", pid, tempCountGreater, pointsCoordinates[tempCountGreater][0],
                //       pointsCoordinates[tempCountGreater][1]);
                tempCountGreater++;
            }
        }/**/
        //printf("(AFTER) ProcessID: %d\tLessEqual = %d\tGreater = %d.\n", processID, tempCountLessEqual, tempCountGreater);

        // Masters create the node and send it to allmaster
        int treeOffset;
        if (pid == 0) {
            // Masters construct the tree
            treeOffset = groups - 1; // This is the index of the start of the nodes on the same height
            //printf("treeOffset = %d\n", treeOffset);
            treeOffset += processID / noProcesses; // Find the exact position of the node on the tree
            //for (int m = 0; m < groups; ++m) {
            tree[treeOffset].radius = median; // Store the radius value of the node
            tree[treeOffset].vp = (float *) malloc(d * sizeof(float));
            for (int i = 0; i < d; ++i) {
                tree[treeOffset].vp[i] = vp[i]; // Store the Vantage Point
            }

            if (processID == 0) {
                // Receive all the nodes of this height
                treeOffset = groups - 1;
                for (int i = 1; i < groups; ++i) {
                    MPI_Recv(&(tree[treeOffset + i].radius), 1, MPI_FLOAT, i * noProcesses, 8, MPI_COMM_WORLD, &mpiStat);
                    tree[treeOffset + i].vp = (float *) malloc(d * sizeof(float));
                    MPI_Recv(tree[treeOffset + i].vp, d, MPI_FLOAT, i * noProcesses, 9, MPI_COMM_WORLD, &mpiStat);
                }
            } else {
                // Send the node to the master (of the world)
                MPI_Send(&median, 1, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
                MPI_Send(vp, d, MPI_FLOAT, 0, 9, MPI_COMM_WORLD);
                free(tree[treeOffset].vp);
            }

            // Free for malloc
            free(countersLE);
            free(countersG);
            free(statuses);
        }

        //Allmaster sends the nodes to all the processes
        treeOffset = groups - 1;
        if (processID == 0) {
            for (int i = 0; i < groups; ++i) { // Iterate all nodes to be sent
                MPI_Bcast(&(tree[treeOffset + i].radius), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                MPI_Bcast(tree[treeOffset + i].vp, d, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //printf("Process %d\tTree height %d,%d\t r=%f\tvp=(%f, %f)\n", processID, loop, i, tree[treeOffset + i].radius,
                //       tree[treeOffset + i].vp[0], tree[treeOffset + i].vp[1]);
            }
        } else {
            for (int i = 0; i < groups; ++i) { // Iterate all nodes to be received
                MPI_Bcast(&(tree[treeOffset + i].radius), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                tree[treeOffset + i].vp = (float *) malloc(d * sizeof(float));
                MPI_Bcast(tree[treeOffset + i].vp, d, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //printf("Process %d\tTree height %d,%d\t r=%f\tvp=(%f, %f)\n", processID, loop, i, tree[treeOffset + i].radius,
                //       tree[treeOffset + i].vp[0], tree[treeOffset + i].vp[1]);
            }
        }

        MPI_Barrier(*communicator);
    }

    // SINGLE PROCESS PART
    MPI_Comm_split(MPI_COMM_WORLD, processID, 0, communicator); // Each process becomes a master with no slaves. Use for mpiFindMedian.
    MPI_Comm_rank(*communicator, &pid);

    int treeOffset = noTotalProcesses - 1;
    //printf("index = %d perProcessSize = %d\n", treeOffset + processID, perProcessSize);
    finishTree(processID, tree, treeOffset + processID, pointsCoordinates, perProcessSize, d, communicator); // Finish the tree locally
    MPI_Barrier(*communicator);

    for (int m = 0; m < noTotalProcesses; ++m) {
        shareTree(processID, m, tree, treeOffset + m, 2 * totalSize - 1, d);
    }

    for (int n = 0; n < totalSize; ++n) {
        int temp = totalSize - 1 + n;
        printf("Process %d tree[%d] -> r = %f\tvp = (%f, %f)\n", processID, temp, tree[temp].radius, tree[temp].vp[0], tree[temp].vp[1]);
    }

    // Free for malloc
    free(distances);
    free(tempDistances);
    for (int k = 0; k < perProcessSize; ++k) {
        free(pointsCoordinates[k]);
        free(lessEqual[k]);
        free(greater[k]);
    }
    free(pointsCoordinates);
    free(lessEqual);
    free(greater);
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Main Median = %f\n",median);

    MPI_Type_free(&mpi_node_type);
    MPI_Finalize();

}
