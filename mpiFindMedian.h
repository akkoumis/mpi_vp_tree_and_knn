//
// Created by Anastasios on 13/12/2018.
//

#ifndef EX2_MPIFINDMEDIAN_H
#define EX2_MPIFINDMEDIAN_H

#include <mpi.h>

void generateNumbers(float *numberPart, int partLength, int processID);

float mpiFindMedian(int processId, int noProcesses, int sizeOfArray, float *distances, MPI_Comm *communicator);

#endif //EX2_MPIFINDMEDIAN_H
