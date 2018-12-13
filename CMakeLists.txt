cmake_minimum_required(VERSION 3.10)
project(ex2)

set(CMAKE_C_STANDARD 11)
find_package(MPI)
include_directories(.)
include_directories(/opt/intel/compilers_and_libraries/linux/mpi/intel64/include)

add_executable(ex2
        mpiFindMedian.c)
target_link_libraries(ex2 ${MPI_C_LIBRARIES})