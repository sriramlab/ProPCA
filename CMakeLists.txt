cmake_minimum_required(VERSION 3.15)
project(ProPCA VERSION 1.0.0 LANGUAGES CXX)
include(GNUInstallDirs)
set(default_build_type "Release")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG=1")


set(CMAKE_CXX_STANDARD 11)
find_package(Threads REQUIRED)

option(PROPCA_ASAN "Enable (debug) address sanitizer build" OFF)
option(PROPCA_TSAN "Enable (debug) Thread sanitizer build" OFF)
add_executable(propca src/propca.cpp src/genotype.cpp src/storage.cpp)
target_include_directories(propca PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(propca Threads::Threads)



if(PROPCA_ASAN)
  set(CMAKE_BUILD_TYPE "Debug")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif()


if(PROPCA_TSAN)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=thread")
endif()


install(TARGETS propca)
