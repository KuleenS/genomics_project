cmake_minimum_required(VERSION 3.22) 
project(computational_genomics)

find_path(RAPIDCSV_INCLUDE_DIRS "rapidcsv.h")

file(GLOB FOUR_RUSSIANS
     "src/four_russians/*.h"
     "src/four_russians/*.cpp"
)

add_executable(alignment src/alignment.cpp src/algos.cpp src/algos.hpp)
add_executable(hirschberg src/hirschberg.cpp)
add_executable(four_russians src/four_russians.cpp ${FOUR_RUSSIANS})

target_include_directories(alignment PRIVATE ${RAPIDCSV_INCLUDE_DIRS})
target_include_directories(hirschberg PRIVATE ${RAPIDCSV_INCLUDE_DIRS})
