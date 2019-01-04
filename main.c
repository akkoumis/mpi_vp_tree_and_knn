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

int d = 2;
int maxK = 4;

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
    float res = sqrtf(sum);
    //printf("Distance to (%f,%f):\t%f\n", end[0], end[1], res);
    return res;
}

void copyPoint(float *new, const float *old, int dimension) {
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
        tree[index].radius = -99;                // Insert node's radius
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
        //printf("tree[%d] r = %f vp = (%f,%f)\n", index, median, tree[index].vp[0], tree[index].vp[1]);
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

void fillTreeWithOtherProcesses(int myProcessID, int otherProcessID, node *tree, int index, int treeSize) {
    if (myProcessID != otherProcessID) {
        tree[index].radius = -otherProcessID - 1;
        tree[index].vp = NULL;
        if ((index * 2 + 2) < treeSize) {
            fillTreeWithOtherProcesses(myProcessID, otherProcessID, tree, index * 2 + 1, treeSize);
            fillTreeWithOtherProcesses(myProcessID, otherProcessID, tree, index * 2 + 2, treeSize);
        }
    }
}

void findLeafsAndCalculateDistance(const node *tree, int startIndex, int *counter, float **leafs, const float *qp, int dimension) {
    if (tree[startIndex].radius == -99) { // It means we found a leaf
        copyPoint(leafs[*counter], tree[startIndex].vp, dimension);
        leafs[*counter][dimension] = calculateDistance(qp, leafs[*counter], dimension); // Add distance from qp to the end
        //printf("Leaf[%d] = (%f, %f) distFromVP = %f\n", *counter, leafs[*counter][0], leafs[*counter][1], leafs[*counter][dimension]);
        (*counter)++;
    } else {
        findLeafsAndCalculateDistance(tree, 2 * startIndex + 1, counter, leafs, qp, dimension);
        findLeafsAndCalculateDistance(tree, 2 * startIndex + 2, counter, leafs, qp, dimension);
    }
}

void getCandidates(const node *tree, int startIndex, int firstStartIndex, int *counterCandidates, float **candidates, const float *qp, const
float kNNDistance) {
    if (tree[startIndex].radius == -99) { // Reached the leafs
        float distanceQPToLeaf = calculateDistance(qp, tree[startIndex].vp, d);
        if (distanceQPToLeaf < kNNDistance) {
            //printf("Found candidate[%d]\n", *counterCandidates);
            if (*counterCandidates > 0) {
                int i = 0;
                //printf("ERROR 1\n");
                //printf("candidate[%d] dist = %f\n", i, calculateDistance(qp, candidates[i], d));
                while (distanceQPToLeaf > calculateDistance(qp, candidates[i], d)) {
                    i++;
                    if (i >= *counterCandidates)
                        break;
                }
                //printf("ERROR 2\n");
                for (int j = *counterCandidates; j > i; j--) {
                    copyPoint(candidates[j], candidates[j - 1], d + 1);
                }
                //printf("Inserting candidate[%d] at index %d\n", *counterCandidates, i);
                copyPoint(candidates[i], tree[startIndex].vp, d);
                candidates[i][d] = distanceQPToLeaf;

            } else {/**/
                //printf("Inserting candidate[%d] at index %d\n", *counterCandidates, *counterCandidates);
                copyPoint(candidates[*counterCandidates], tree[startIndex].vp, d);
                candidates[*counterCandidates][d] = distanceQPToLeaf;
            }
            (*counterCandidates)++;
            // TODO Maybe change knnDistance, by removing const and converting to pointer
        }
    } else {
        // TODO Prune the tree accordingly
        if ((2 * startIndex + 1) != firstStartIndex && tree[startIndex].radius >= (calculateDistance(tree[startIndex].vp, qp, d) - kNNDistance))
            getCandidates(tree, 2 * startIndex + 1, firstStartIndex, counterCandidates, candidates, qp, kNNDistance);
        if ((2 * startIndex + 2) != firstStartIndex && tree[startIndex].radius <= (calculateDistance(tree[startIndex].vp, qp, d) + kNNDistance))
            getCandidates(tree, 2 * startIndex + 2, firstStartIndex, counterCandidates, candidates, qp, kNNDistance);
    }
}

void getCandidatesOthers(const node *tree, int startIndex, int *counterCandidates, float **candidates, const float *qp, const float kNNDistance) {
    //printf("Searching node %d\n", startIndex);
    if (tree[startIndex].radius == -99) { // Reached the leafs
        float distanceQPToLeaf = calculateDistance(qp, tree[startIndex].vp, d);
        if (distanceQPToLeaf < kNNDistance) {
            //printf("Found candidate[%d]\n", *counterCandidates);
            if (*counterCandidates > 0) {
                int i = 0;
                while (distanceQPToLeaf > calculateDistance(qp, candidates[i], d)) {
                    i++;
                    if (i >= *counterCandidates)
                        break;
                }
                for (int j = *counterCandidates; j > i; j--) {
                    copyPoint(candidates[j], candidates[j - 1], d + 1);
                }
                //printf("Inserting candidate[%d] at index %d\n", *counterCandidates, i);
                copyPoint(candidates[i], tree[startIndex].vp, d);
                candidates[i][d] = distanceQPToLeaf;
            } else {/**/
                //printf("Inserting candidate[%d] at index %d\n", *counterCandidates, *counterCandidates);
                copyPoint(candidates[*counterCandidates], tree[startIndex].vp, d);
                candidates[*counterCandidates][d] = distanceQPToLeaf;
            }
            (*counterCandidates)++;
            // TODO Maybe change knnDistance, by removing const and converting to pointer
        }
    } else {
        // TODO Prune the tree accordingly
        if (tree[startIndex].radius >= (calculateDistance(tree[startIndex].vp, qp, d) - kNNDistance))
            getCandidatesOthers(tree, 2 * startIndex + 1, counterCandidates, candidates, qp, kNNDistance);
        if (tree[startIndex].radius <= (calculateDistance(tree[startIndex].vp, qp, d) + kNNDistance))
            getCandidatesOthers(tree, 2 * startIndex + 2, counterCandidates, candidates, qp, kNNDistance);
    }
}


int compareBasedOnDistance(const void *a, const void *b) {
    float x = (*(const float **) a)[d];
    float y = (*(const float **) b)[d];
    int temp;
    temp = x < y ? -1 : x == y ? 0 : 1;
    //printf("Compare = %d\n", temp);
    return temp;
}

int compareBasedOnCandidatesDistance(const void *a, const void *b) {
    float x = (*(const float **) a)[d];
    float y = (*(const float **) b)[d];
    int temp;
    temp = x < y ? -1 : x == y ? 0 : 1;
    //printf("Compare = %d\n", temp);
    return temp;
}

void findClosestFrontierDistance(const node *tree, const int startIndex, const int qpIndex, float *closestFrontierDistance) {
    float distanceToFrontier;
    int parentIndex = (startIndex - 1) / 2;
    if (parentIndex >= 0) {
        if (startIndex % 2 == 1) {// It means the node tree[startIndex] is a left child => inside of the radius
            distanceToFrontier = tree[parentIndex].radius - calculateDistance(tree[qpIndex].vp, tree[parentIndex].vp, d);
        } else {
            distanceToFrontier = calculateDistance(tree[qpIndex].vp, tree[parentIndex].vp, d) - tree[parentIndex].radius;
        }
        //printf("Pid %pid distanceFrontierOfIndex[%d] = %f\n",pid, startIndex, distanceToFrontier);
        if (distanceToFrontier < *closestFrontierDistance)
            *closestFrontierDistance = distanceToFrontier;
        if (parentIndex > 0)
            findClosestFrontierDistance(tree, parentIndex, qpIndex, closestFrontierDistance);
    }
}


int main(int argc, char **argv) {
    MPI_Status mpiStat;


    int processID, pid, noProcesses, noTotalProcesses;
    int totalSize, perProcessSize, loop, master, groups, perGroupSize; // Size = # of elems

    float *distances, **pointsCoordinates, **lessEqual, **greater;;
    totalSize = 128;

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
    node *tree2 = (node *) malloc((2 * totalSize - 1) * sizeof(node));
//printf("tree size = %d\n", (2 * totalSize - 1));

    // Each process read simultaneously its data, from the file
    FILE *fp;
    char buff[255], buff2[255];

    long int offset = ((long int) perProcessSize) * 2 * processID; // *2 -> because x and y are on different lines
    fp = fopen("cities.txt", "r");
    for (long int i = 0; i < offset; ++i) {
        fscanf(fp, "%s", buff); // Moving stream according to offset
    }

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

    char *endptr;
    for (int i = 0; i < perProcessSize; ++i) {
        for (int j = 0; j < d; ++j) {
            fscanf(fp, "%s", buff2);
            pointsCoordinates[i][j] = strtof(buff2, &endptr);
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
        float median = mpiFindMedian(processID, noProcesses, perGroupSize, tempDistances, communicator);

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
        //MPI_Barrier(*communicator);

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

    int treeOffset = noTotalProcesses - 1; // Points to the start of the local subtrees
    //printf("index = %d perProcessSize = %d\n", treeOffset + processID, perProcessSize);
    finishTree(processID, tree, treeOffset + processID, pointsCoordinates, perProcessSize, d, communicator); // Finish the tree locally
    MPI_Barrier(MPI_COMM_WORLD);

    for (int m = 0; m < 2 * totalSize - 1; ++m) {
        tree2[m] = tree[m];
        //printf("PID %d tree[%d] -> r = %f\n", processID, m, tree[m].radius);
    }

    for (int m = 0; m < noTotalProcesses; ++m) {
        shareTree(processID, m, tree2, treeOffset + m, 2 * totalSize - 1, d); // Share the tree, for debugging
        fillTreeWithOtherProcesses(processID, m, tree, treeOffset + m, 2 * totalSize - 1); // Fill the blanks with holder's processID
    }




    /*for (int n = 0; n < 2 * totalSize - 1; ++n) {
        int temp = totalSize - 1 + n;
        if (tree[n].radius >= 0 || tree[n].radius == -99)//(tree[n].vp != NULL)
            printf("Process %d tree[%d] -> r = %f\tvp = (%f, %f)\n", processID, n, tree[n].radius, tree[n].vp[0], tree[n].vp[1]);
        else
            printf("Process %d tree[%d] -> r = %f\n", processID, n, tree[n].radius);
    }*/

    //MPI_Barrier(MPI_COMM_WORLD);

    // k-ΝΝ
    float ***neighbors = (float ***) malloc(perProcessSize * sizeof(float **)); // To store the neighbors of each local point
    for (int l = 0; l < perProcessSize; ++l) { // perProcessSize array needed for each point to store its kNN
        neighbors[l] = (float **) malloc(maxK * sizeof(float *)); // k sized array needed for the neighbors
        for (int j = 0; j < maxK; ++j) {
            neighbors[l][j] = (float *) malloc((d + 1) * sizeof(float));
        }
    }

    float distanceToFrontier;
    //for (int i = 1; i <= 1; ++i) { // k = 2^[1:8]
    int i = 2;
    int k = (int) pow(2, i); // Calculate the k

    maxK = k;

    MPI_Barrier(MPI_COMM_WORLD);
    for (int searchProcess = 0; searchProcess < noTotalProcesses; ++searchProcess) { // TODO process < noTotalProcesses
        MPI_Barrier(MPI_COMM_WORLD);
        if (searchProcess == processID) { // It means that this process execute the search
            int startIndex, parentIndex, qpIndex;

            treeOffset = totalSize - 1 + processID * perProcessSize; // Points to the start of the leafs.
            printf("\nProcess %d is searching at [%d]\n", processID, treeOffset);

            //float ***pointsToBeSent = (float ***) malloc(noTotalProcesses * sizeof(float **));
            int pointsToBeSentIndices[perProcessSize][noTotalProcesses];
            int *countersPointsToBeSent = (int *) malloc(noTotalProcesses * sizeof(int));
            for (int n = 0; n < noTotalProcesses; ++n) {
                countersPointsToBeSent[n] = 0; // How many points to send to each process
            }

            for (int j = 0; j < 1; ++j) { // TODO j<perProcessSize
                qpIndex = treeOffset + j;
                for (int n = 0; n < noTotalProcesses; ++n) {
                    pointsToBeSentIndices[j][n] = 0; // Is the point j to be sent to process n?
                }

                startIndex = treeOffset + j;
                for (int l = 0; l < i; ++l) {
                    startIndex = (startIndex - 1) / 2; // Find the start of the subtree that has as many
                }
                //printf("startIndex = %d\n",startIndex);
                parentIndex = (startIndex - 1) / 2;

                int counterLeaf[] = {0};
                findLeafsAndCalculateDistance(tree, startIndex, counterLeaf, neighbors[j], tree[qpIndex].vp, d);
                qsort(neighbors[j], k, sizeof(*neighbors), compareBasedOnDistance);
                /*for (int m = 0; m < k; ++m) {
                    printf("QNN[%d] = (%f, %f) distFromVP = %f\n", m, neighbors[j][m][0], neighbors[j][m][1], neighbors[j][m][d]);
                }*/

                if (startIndex % 2 == 1) {// It means the node tree[startIndex] is a left child => inside of the radius
                    distanceToFrontier = tree[parentIndex].radius - calculateDistance(tree[qpIndex].vp, tree[parentIndex].vp, d);
                } else {
                    distanceToFrontier = calculateDistance(tree[qpIndex].vp, tree[parentIndex].vp, d) - tree[parentIndex].radius;
                }

                findClosestFrontierDistance(tree, startIndex, qpIndex, &distanceToFrontier);

                printf("PID = %d\tkNNDist = %f\tfrontierDist = %f\n", processID, neighbors[j][k - 1][d], distanceToFrontier);
                if (neighbors[j][k - 1][d] > distanceToFrontier) { // If furthest kNN is further than the frontier
                    printf("kΝΝ pid=%d is further than the frontier\n", processID);

                    // At first we check the local subtree (of this process)
                    int otherProcessesNeeded = 0;
                    float distanceToLocalFrontier;
                    int localTreeStartIndex = noTotalProcesses - 1 + processID;
                    int parentLocal = (localTreeStartIndex - 1) / 2;


                    float **candidatePoints = (float **) malloc(perProcessSize * sizeof(float *));
                    for (int l = 0; l < perProcessSize; ++l) {
                        candidatePoints[l] = (float *) malloc((d + 1) * sizeof(float));
                    }
                    int counterCandidates[] = {0};
                    //printf("At PID = %d\tlocalTreeStartIndex = %d\tstartIndex = %d\n", processID, localTreeStartIndex, startIndex);
                    if (localTreeStartIndex != startIndex) {
                        printf("At PID = %d searching for candidates starting at node %d\n", processID, localTreeStartIndex);
                        //printf("ERROR 1\n");

                        getCandidates(tree, localTreeStartIndex, startIndex, counterCandidates, candidatePoints, tree[qpIndex].vp,
                                      neighbors[j][k - 1][d]);
                        //printf("ERROR 2\n");
                    }
                    for (int m = 0; m < *counterCandidates; ++m) {
                        //printf("Candidates[%d] @ pid %d = (%f, %f) distFromVP = %f\n", m, processID, candidatePoints[m][0],
                        //       candidatePoints[m][1], candidatePoints[m][d]);
                    }

                    /*for (int n = 0; n < *counterCandidates - 1; ++n) {
                        for (int l = n + 1; l < *counterCandidates; ++l) {
                            if (candidatePoints[n][d] > candidatePoints[l][d]) {
                                float temp[d + 1];
                                copyPoint(temp, candidatePoints[n], d + 1);
                                copyPoint(candidatePoints[n], candidatePoints[l], d + 1);
                                copyPoint(candidatePoints[l], temp, d + 1);
                            }
                        }
                    }
                    qsort(candidatePoints, *counterCandidates, sizeof(*neighbors), compareBasedOnCandidatesDistance);
                    for (int m = 0; m < *counterCandidates; ++m) {
                        printf("QCandidates[%d] @ pid %d = (%f, %f) distFromVP = %f\n", m, processID, candidatePoints[m][0],
                               candidatePoints[m][1], candidatePoints[m][d]);
                    }*/


                    // TODO integrate candidate points
                    int indexNeighbors, indexCandidates = 0;
                    for (indexNeighbors = 0; indexNeighbors < k; ++indexNeighbors) {
                        if (neighbors[j][indexNeighbors][d] > candidatePoints[indexCandidates][d] && indexCandidates < *counterCandidates) {
                            for (int l = k - 1; l > indexNeighbors; l--) {
                                copyPoint(neighbors[j][l], neighbors[j][l - 1], d + 1);
                            }
                            float temp[d + 1];
                            copyPoint(temp, neighbors[j][indexNeighbors], d + 1);
                            copyPoint(neighbors[j][indexNeighbors], candidatePoints[indexCandidates], d + 1);
                            copyPoint(candidatePoints[indexCandidates], temp, d + 1);
                            indexCandidates++;
                        }
                    }


                    if (localTreeStartIndex % 2 == 1) {// It means the node of the local subtree is a left child => inside of the radius
                        distanceToLocalFrontier =
                                tree[parentLocal].radius - calculateDistance(tree[qpIndex].vp, tree[parentLocal].vp, d);
                    } else {
                        distanceToLocalFrontier =
                                calculateDistance(tree[qpIndex].vp, tree[parentLocal].vp, d) - tree[parentLocal].radius;
                    }

                    findClosestFrontierDistance(tree, localTreeStartIndex, qpIndex, &distanceToLocalFrontier);

                    if (neighbors[j][k - 1][d] > distanceToLocalFrontier) {
                        otherProcessesNeeded = 1;
                    }

                    // TODO fill this with communications between processes
                    if (otherProcessesNeeded) {
                        printf("Other processes are needed at pid = %d\n", processID);

                        for (int l = 0; l < noTotalProcesses; ++l) { //TODO nototalpro
                            if (l != processID) {
                                int count = 1;
                                MPI_Send(&count, 1, MPI_INT, l, 10, MPI_COMM_WORLD); // Send signal
                                MPI_Send(tree[qpIndex].vp, d, MPI_FLOAT, l, 11, MPI_COMM_WORLD); // Send qp
                                MPI_Send(&neighbors[j][k - 1][d], 1, MPI_FLOAT, l, 12, MPI_COMM_WORLD); // Send kNN distance
                                float temp[d + 1];
                                MPI_Recv(&count, 1, MPI_INT, l, 13, MPI_COMM_WORLD, &mpiStat); // Receive candidates counter
                                int recvdSize = count * (d + 1);
                                float *recvdCandidates = (float *) malloc(recvdSize * sizeof(float));
                                MPI_Recv(recvdCandidates, recvdSize, MPI_FLOAT, l, 14, MPI_COMM_WORLD, &mpiStat); // Receive candidates
                                for (int m = 0; m < count; ++m) {
                                    copyPoint(candidatePoints[m], &recvdCandidates[m * (d + 1)], d + 1);
                                }/**/


                                int startOffset = noTotalProcesses - 1;
                                //printf("Other searching PID %d at index = %d for kNNDistance = %f\n", processID, startOffset + processID,
                                //        neighbors[j][k - 1][d]);

                                int oCounterCandidates[] = {0};
                                float **oCandidatePoints = (float **) malloc(perProcessSize * sizeof(float *));
                                for (int p = 0; p < perProcessSize; ++p) {
                                    oCandidatePoints[p] = (float *) malloc((d + 1) * sizeof(float));
                                }
                                /*getCandidatesOthers(tree, startOffset + l, oCounterCandidates, oCandidatePoints, tree[qpIndex].vp,
                                                    neighbors[j][k - 1][d]);
                                printf("Other searching PID %d at index = %d. Found %d cands\n", processID, startOffset + processID,
                                       *oCounterCandidates);
                                for (int m = 0; m < count; ++m) {
                                    printf("PID %d from %d candidate[%d] = (%f, %f) distFromVP = %f\n", processID, l, m, oCandidatePoints[m][0],
                                           oCandidatePoints[m][1], oCandidatePoints[m][d]);
                                }*/


                                int oIndexNeighbors, oIndexCandidates = 0;
                                for (oIndexNeighbors = 0; oIndexNeighbors < k; ++oIndexNeighbors) {
                                    if (neighbors[j][oIndexNeighbors][d] > candidatePoints[oIndexCandidates][d] && oIndexCandidates < count) {
                                        for (int oL = k - 1; oL > oIndexNeighbors; oL--) {
                                            copyPoint(neighbors[j][oL], neighbors[j][oL - 1], d + 1);
                                        }
                                        float otemp[d + 1];
                                        copyPoint(otemp, neighbors[j][oIndexNeighbors], d + 1);
                                        copyPoint(neighbors[j][oIndexNeighbors], candidatePoints[oIndexCandidates], d + 1);
                                        copyPoint(candidatePoints[oIndexCandidates], otemp, d + 1);
                                        oIndexCandidates++;
                                    }
                                }

                                free(recvdCandidates);
                                for (int p = 0; p < perProcessSize; ++p) {
                                    free(oCandidatePoints[p]);
                                }
                                free(oCandidatePoints);
                            }
                        }

                    }



                    //free(pointsToBeSent);
                    free(countersPointsToBeSent);
                    for (int l = 0; l < perProcessSize; ++l) {
                        free(candidatePoints[l]);
                    }
                    free(candidatePoints);
                }
            }

            int finishCode = -1;
            for (int n = 0; n < noTotalProcesses; ++n) {
                if (n != searchProcess) {
                    //printf("%d is sending finish signal to %d\n", processID, n);
                    MPI_Send(&finishCode, 1, MPI_INT, n, 10, MPI_COMM_WORLD);
                }
            }
//MPI_Barrier(MPI_COMM_WORLD);
        } else { // It means this process helps the one that searches
            int temp = 0; // To be renamed to counter
            MPI_Recv(&temp, 1, MPI_INT, searchProcess, 10, MPI_COMM_WORLD, &mpiStat);
            //printf("%d received finish signal from %d\n", processID, searchProcess);
            while (temp != -1) {
                // TODO search for candidates
                float qp[d + 1];
                float kNNDistance;
                MPI_Recv(qp, d, MPI_FLOAT, searchProcess, 11, MPI_COMM_WORLD, &mpiStat); // Receive qp
                MPI_Recv(&kNNDistance, 1, MPI_FLOAT, searchProcess, 12, MPI_COMM_WORLD, &mpiStat); // Receive kNNDistance
                float **candidatePoints = (float **) malloc(perProcessSize * sizeof(float *));
                for (int l = 0; l < perProcessSize; ++l) {
                    candidatePoints[l] = (float *) malloc((d + 1) * sizeof(float));
                }
                int counterCandidates[] = {0};

                int startOffset = noTotalProcesses - 1;

                getCandidates(tree, startOffset + processID, startOffset + searchProcess, counterCandidates, candidatePoints, qp, kNNDistance);

                MPI_Send(counterCandidates, 1, MPI_INT, searchProcess, 13, MPI_COMM_WORLD); // Send candidates counter
                int sentSize = ((d + 1) * (*counterCandidates));
                float *sentCandidates = (float *) malloc(sentSize * sizeof(float));
                for (int j = 0; j < *counterCandidates; ++j) {
                    copyPoint(&sentCandidates[j * (d + 1)], candidatePoints[j], d + 1);
                }
                MPI_Send(sentCandidates, sentSize, MPI_FLOAT, searchProcess, 14, MPI_COMM_WORLD); // Send candidates


                for (int l = 0; l < perProcessSize; ++l) {
                    free(candidatePoints[l]);
                }
                free(candidatePoints);
                free(sentCandidates);
                MPI_Recv(&temp, 1, MPI_INT, searchProcess, 10, MPI_COMM_WORLD, &mpiStat);/**/
            }
            //MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    //}
    MPI_Barrier(MPI_COMM_WORLD); // For timing



    // VALIDATION

    /*for (int m = 0; m < noTotalProcesses; ++m) {
        shareTree(processID, m, tree, treeOffset + m, 2 * totalSize - 1, d); // Share the tree, for validation
    }


    for (int n = 0; n < 2 * totalSize - 1; ++n) {
        int temp = totalSize - 1 + n;
        if (tree[n].radius >= 0 || tree[n].radius == -99)//(tree[n].vp != NULL)
            printf("Process %d sharedTree[%d] -> r = %f\tvp = (%f, %f)\n", processID, n, tree[n].radius, tree[n].vp[0], tree[n].vp[1]);
        else
            printf("Process %d sharedTree[%d] -> r = %f\n", processID, n, tree[n].radius);
    }*/

    float ***validNeighbors = (float ***) malloc(perProcessSize * sizeof(float **)); // To store the neighbors of each local point
    for (int l = 0; l < perProcessSize; ++l) { // perProcessSize array needed for each point to store its kNN
        validNeighbors[l] = (float **) malloc(totalSize * sizeof(float *)); // k sized array needed for the neighbors
        for (int j = 0; j < totalSize; ++j) {
            validNeighbors[l][j] = (float *) malloc((d + 1) * sizeof(float));
        }
    }

    for (int j = 0; j < 1; ++j) { // TODO j<perProcessSize for every point of the process
        int o = processID * perProcessSize + j;
        int counterLeaf[] = {0};
        //printf("ERROR 1\n");
        findLeafsAndCalculateDistance(tree2, 0, counterLeaf, validNeighbors[j], tree[totalSize - 1 + o].vp, d);
        qsort(validNeighbors[j], totalSize, sizeof(*validNeighbors), compareBasedOnDistance);
        int valid = 1;
        for (int i = 0; i < maxK; ++i) {
            printf("PID %d (%f,%f) %f\t(%f,%f) %f\n", processID, neighbors[j][i][0], neighbors[j][i][1], neighbors[j][i][d],
                   validNeighbors[j][i][0], validNeighbors[j][i][1], validNeighbors[j][i][d]);
            for (int k = 0; k < d + 1; ++k) {
                if (validNeighbors[j][i][k] != neighbors[j][i][k])
                    valid = 0;
            }
        }
        if (valid)
            printf("NN PID %d PASSED\n", processID);
        else
            printf("NN PID %d VALIDATION FAILED\n", processID);
    }

    // Free neighbor array
    for (int l = 0; l < perProcessSize; ++l) { // perProcessSize array needed for each point to store its kNN
        for (int j = 0; j < maxK; ++j) {
            free(neighbors[l][j]);
        }
        free(neighbors[l]);
    }
    free(neighbors);
    free(tree);
    free(tree2);


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

    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Process %d EXITING...\n", processID);
    MPI_Type_free(&mpi_node_type);
    MPI_Finalize();

}
