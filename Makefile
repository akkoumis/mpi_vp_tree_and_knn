ICC = mpiicc
ICCFLAGS = -I /usr/include/x86_64-linux-gnu/c++/8/
#OPENMPICCFLAGS = -qopenmp
#PTHREADSICCFLAGS = -pthread

#GCC = gcc
#GCCFLAGS = -fcilkplus -O3
#OPENMPGCCCFLAGS = -fopenmp

CMAIN = mpiFindMedian
SOURCE = mpiFindMedian
NPROCS = 4


all: mpiFindMedian.o
	$(ICC) $(ICCFLAGS) $(OPENMPICCFLAGS) $^ -o $(CMAIN)

mpiFindMedian.o: mpiFindMedian.c
	$(ICC) -c $(ICCFLAGS) $^
#qsort-parallel-cilk.o: qsort-parallel-cilk.c
#	$(ICC) -c $(ICCFLAGS) $^
#qsort-parallel-openmp.o: qsort-parallel-openmp.c
#	$(ICC) -c $(ICCFLAGS) $(OPENMPICCFLAGS) $^
#qsort-parallel-pthreads.o: qsort-parallel-pthreads.c
#	$(ICC) -c $(ICCFLAGS) $(PTHREADICCFLAGS) $^
#%.o: %.c
#	$(ICC) -c $(ICCFLAGS) $^
run:	$(SOURCE).o
	# mpiexec -machinefile mpd.hosts -n $(NPROCS) ./$(SOURCE).out
	mpirun -n $(NPROCS) ./$(SOURCE) 512

clean:
	rm -f *.o *~ $(CMAIN)
