cmake_minimum_required(VERSION 3.27)
project(listupdate)

set(CMAKE_CXX_STANDARD 23)

add_executable(listupdate main.cpp
        memory.hpp
        memory_test.cpp
        permutations.hpp
        permutation_test.cpp
        canonical_ordering.hpp)
