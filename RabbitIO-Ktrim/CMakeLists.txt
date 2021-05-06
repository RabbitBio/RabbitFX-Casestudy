cmake_minimum_required(VERSION 3.0)

project(RabbitIO-Ktrim)
#openmp
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERB")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH .)
endif()
#link

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -ffast-math ${CMAKE_CXX_FLAGS}")

add_subdirectory(src/io)
#message(${SOURCE_LIST})
#include_directories(./src)

add_executable (${PROJECT_NAME} src/main.cpp)
#add_executable (RabbitIO ${SOURCE_LIST})

target_link_libraries(${PROJECT_NAME} z)
target_link_libraries(${PROJECT_NAME} io_lib)
