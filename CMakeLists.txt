cmake_minimum_required(VERSION 3.13)
project(Learned_Rank_Select_ALENEX21)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread -fopenmp")
include_directories(include lib/la_vector/include lib/ds2i)

add_subdirectory(lib/ds2i)
add_subdirectory(lib/la_vector)
add_subdirectory(tests)


add_executable(Comparison comparison.cpp)
target_link_libraries(Comparison sdsl)

add_executable(Generate_datasets generate_datasets.cpp)
target_link_libraries(Generate_datasets sdsl)
