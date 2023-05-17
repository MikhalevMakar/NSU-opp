cmake_minimum_required(VERSION 3.16.3)
project(Task5)

set(CMAKE_CXX_STANDARD 14)

add_executable(Task5 main.cpp ConcurrentQueue.h)

find_package(MPI REQUIRED)

include_directories("${MPI_CXX_INCLUDE_DIRS}")

target_link_libraries(Task5 ${MPI_LIBRARIES})
