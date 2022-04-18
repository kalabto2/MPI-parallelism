cmake_minimum_required(VERSION 3.21)
project(process_paralellism_new)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fopenmp")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Add MPI Package to Project
find_package(MPI REQUIRED)

add_executable(process_paralellism_new main.cpp)

# Add libraries for code completion and compiling
target_link_libraries(process_paralellism_new PUBLIC MPI::MPI_CXX)
