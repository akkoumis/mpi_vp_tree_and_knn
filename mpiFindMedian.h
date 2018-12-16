//
// Created by Anastasios on 13/12/2018.
//

#ifndef EX2_MPIFINDMEDIAN_H
#define EX2_MPIFINDMEDIAN_H

#include <mpi.h>

void generateNumbers(float *numberPart, int partLength, int processID);

void mpiFindMedian(int processId, int master, int noProcesses, int sizeOfArray, float *distances, int loop,
                   MPI_Comm *communicator);

#endif //EX2_MPIFINDMEDIAN_H
