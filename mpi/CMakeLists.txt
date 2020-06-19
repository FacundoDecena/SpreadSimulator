cmake_minimum_required(VERSION 3.16)
project(mpi C)

# added optimization
set(CMAKE_C_FLAGS "-O3")
# added -fopenmp
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_C_STANDARD 99)

add_executable(mpi main.c)