cmake_minimum_required(VERSION 3.16.3)
project(Version1)

set(CMAKE_CXX_STANDARD 14)

add_executable(Version1 main.cpp)

find_package(MPI REQUIRED)

include_directories("${MPI_CXX_INCLUDE_DIRS}")

target_link_libraries(Version1 ${MPI_LIBRARIES})
