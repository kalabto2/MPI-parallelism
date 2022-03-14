cmake_minimum_required(VERSION 3.21)
project(task_parallelism)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")

include_directories(src)

add_executable(task_parallelism
        src/Graph.cpp
        src/Graph.h
        src/Instance_handler.cpp
        src/Instance_handler.h
        src/main.cpp)

