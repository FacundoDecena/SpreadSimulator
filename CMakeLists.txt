cmake_minimum_required(VERSION 3.16)
project(COVID_19 C)

# added -fopenmp
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

set(CMAKE_C_STANDARD 99)

add_executable(COVID_19 main.c)