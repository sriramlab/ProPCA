cmake_minimum_required(VERSION 3.15)
project(ProPCA VERSION 1.0.0 LANGUAGES CXX)

set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 11)
find_package(Threads REQUIRED)

add_executable(propca src/propca.cpp src/genotype.cpp src/storage.cpp)
target_include_directories(propca PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(propca Threads::Threads)
