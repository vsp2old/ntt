cmake_minimum_required (VERSION 3.0)

project (NTT)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
add_definitions("-g -O2 -Wall -std=c++11")

add_executable(ntt ntt.cpp)
add_executable(ntt1 ntt1.cpp gmath.cpp)
add_executable(ntt2 ntt2.cpp gmath.cpp)
add_executable(ntt3 ntt3.cpp gmath.cpp)
