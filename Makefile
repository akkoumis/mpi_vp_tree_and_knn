ICC = mpiicc
ICCFLAGS = -I /usr/include/x86_64-linux-gnu/c++/8/
#OPENMPICCFLAGS = -qopenmp
#PTHREADSICCFLAGS = -pthread

#GCC = gcc
#GCCFLAGS = -fcilkplus -O3
#OPENMPGCCCFLAGS = -fopenmp

CMAIN = ex2
NPROCS = 1


all: main.o mpiFindMedian.o
	$(ICC) $(ICCFLAGS) $^ -o $(CMAIN)

main.o: main.c
	$(ICC) -c $(ICCFLAGS) $^
mpiFindMedian.o: mpiFindMedian.c
	$(ICC) -c $(ICCFLAGS) $^
#qsort-parallel-openmp.o: qsort-parallel-openmp.c
#	$(ICC) -c $(ICCFLAGS) $(OPENMPICCFLAGS) $^
#qsort-parallel-pthreads.o: qsort-parallel-pthreads.c
#	$(ICC) -c $(ICCFLAGS) $(PTHREADICCFLAGS) $^
#%.o: %.c
#	$(ICC) -c $(ICCFLAGS) $^
run:	#$(CMAIN).o
	# mpiexec -machinefile mpd.hosts -n $(NPROCS) ./$(SOURCE).out
	mpirun -n $(NPROCS) ./$(CMAIN) 512

clean:
	rm -f *.o *~ $(CMAIN)
