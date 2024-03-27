cmake_minimum_required(VERSION 3.27)
project(listupdate)

set(CMAKE_CXX_STANDARD 23)

add_executable(listupdate
        main.cpp
        memory_pairs.hpp
        permutations.hpp
        common.hpp
        graph.hpp
        algorithm.hpp
        bellman_ford.hpp
        memory_bitfield.hpp
        implicit_graph.hpp
        implicit_bellman_ford.hpp)
add_executable(permutation_test permutation_test.cpp)
