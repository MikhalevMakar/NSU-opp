cmake_minimum_required(VERSION 3.16.3)
project(Task4)

set(CMAKE_CXX_STANDARD 14)

add_executable(Task4 main.cpp)

find_package(MPI REQUIRED)

include_directories("${MPI_CXX_INCLUDE_DIRS}")

target_link_libraries(Task4 ${MPI_LIBRARIES})