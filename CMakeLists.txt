cmake_minimum_required(VERSION 3.15)
project(ProPCA VERSION 1.0.0 LANGUAGES CXX)
include(GNUInstallDirs)
set(default_build_type "Release")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG=1")

cmake_policy(SET CMP0069 NEW)
include(CheckIPOSupported)
check_ipo_supported(RESULT iporesult)

set(CMAKE_CXX_STANDARD 11)
find_package(Threads REQUIRED)
set(ARCH "haswell" CACHE STRING "Architecture to tell gcc to optimize for (-march)")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=${ARCH} -mtune=${ARCH}")
option(PROPCA_ASAN "Enable (debug) address sanitizer build" OFF)
option(PROPCA_TSAN "Enable (debug) Thread sanitizer build" OFF)
add_executable(propca src/propca.cpp src/genotype.cpp src/storage.cpp)
target_include_directories(propca PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(propca Threads::Threads)
if(iporesult)
  message (STATUS "IPO supported for project ${PROJECT_NAME}")
  set_property(TARGET propca PROPERTY CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")	
endif()



if(PROPCA_ASAN)
  set(CMAKE_BUILD_TYPE "Debug")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif()


if(PROPCA_TSAN)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=thread")
endif()


install(TARGETS propca)
