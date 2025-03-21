cmake_minimum_required(VERSION 3.27)
project(listupdate)

set(CMAKE_CXX_STANDARD 23)


add_executable(algorithm-list-size-4 wf/algorithm_against_opt.cpp)

target_compile_definitions(algorithm-list-size-4 PUBLIC COMP_RATIO=3.0)
target_compile_definitions(algorithm-list-size-4 PUBLIC TSIZE=4)

add_executable(algorithm-list-size-5 wf/algorithm_against_opt.cpp)

target_compile_definitions(algorithm-list-size-5 PUBLIC COMP_RATIO=3.0)
target_compile_definitions(algorithm-list-size-5 PUBLIC TSIZE=5)

add_executable(algorithm-list-size-6 wf/algorithm_against_pairwise_opt.cpp)

target_compile_definitions(algorithm-list-size-6 PUBLIC COMP_RATIO=3.0)
target_compile_definitions(algorithm-list-size-6 PUBLIC TSIZE=6)


add_executable(lowerbound-work-function-class wf/workfunction.cpp)
target_compile_definitions(lowerbound-work-function-class PUBLIC COMP_RATIO=3.1)
target_compile_definitions(lowerbound-work-function-class PUBLIC TSIZE=5)

add_executable(workfunction-reach wf/workfunction_decide.cpp)
add_executable(workfunction_verify wf/workfunction_verify.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(algorithm-list-size-4 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(algorithm-list-size-5 PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(algorithm-list-size-6 PUBLIC OpenMP::OpenMP_CXX)
endif()
