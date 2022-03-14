cmake_minimum_required(VERSION 3.21)
project(sequential)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

add_executable(sequential
        src/Graph.cpp
        src/Graph.h
        src/Instance_handler.cpp
        src/Instance_handler.h
        src/main.cpp)
