cmake_minimum_required(VERSION 3.22) 
project(computational_genomics)

find_path(RAPIDCSV_INCLUDE_DIRS "rapidcsv.h")

add_executable(alignment src/alignment.cpp src/algos.cpp src/algos.hpp)
add_executable(hirschberg src/hirschberg.cpp)

# add_executable(block_alignment src/block_based_alignment.cpp src/algos.cpp src/algos.hpp)


target_include_directories(alignment PRIVATE ${RAPIDCSV_INCLUDE_DIRS})
target_include_directories(hirschberg PRIVATE ${RAPIDCSV_INCLUDE_DIRS})

# target_include_directories(block_alignment PRIVATE ${RAPIDCSV_INCLUDE_DIRS})

