cmake_minimum_required(VERSION 2.8.12)
project(ProPCA)

include_directories(include/)

set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lrt")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})

add_executable(propca src/propca.cpp src/genotype.cpp src/storage.cpp)
