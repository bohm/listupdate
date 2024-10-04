cmake_minimum_required(VERSION 3.27)
project(listupdate)

set(CMAKE_CXX_STANDARD 23)

add_executable(listupdate
        main.cpp
        memory_pairs.hpp
        old_perm_functions.hpp
        common.hpp
        graph.hpp
        algorithm.hpp
        bellman_ford.hpp
        memory_bitfield.hpp
        implicit_graph.hpp
        implicit_bellman_ford.hpp
        storage.hpp
        wf_manager.hpp
        memory_perm.hpp
        iteration_over_memory.hpp
        memory_cartesian.hpp
        graphs/graph_bipartite_mru.hpp
)
add_executable(permutation_test permutation_test.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(listupdate PUBLIC OpenMP::OpenMP_CXX)
endif()
add_executable(workfunction wf/workfunction.cpp
        wf/game_graph.hpp
        wf/double_zobrist.hpp)
add_executable(reachablewf wf/reachablewf.cpp
        wf/double_zobrist.hpp
        data_structures/char_flat_set.hpp)
add_executable(graph_test graph_test.cpp
        graphs/bip_mru_potentials.hpp)
