cmake_minimum_required(VERSION 3.16.3)
project(Task3)

set(CMAKE_CXX_STANDARD 14)

add_executable(Task3 main.cpp)

find_package(MPI REQUIRED)

include_directories("${MPI_CXX_INCLUDE_DIRS}")

target_link_libraries(Task3 ${MPI_LIBRARIES})