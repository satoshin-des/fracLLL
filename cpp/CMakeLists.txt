cmake_minimum_required(VERSION 3.12)
project(fracLLL LANGUAGES CXX)

# set compiler to g++
if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER g++)
endif()

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add compile-options
add_compile_options(
    -O3
    -march=native
    -funroll-loops
    -fomit-frame-pointer
    -fstrict-aliasing
    -ffast-math
    -flto
)

# include directory
include_directories(include)

# collect source files for library
file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")

# build shared library
add_library(fracLLL SHARED ${LIB_SOURCES})

# find NTL
find_library(NTL_LIBRARY NAMES ntl)
if(NOT NTL_LIBRARY)
    message(FATAL_ERROR "NTL library not found")
endif()
target_link_libraries(fracLLL PRIVATE ${NTL_LIBRARY})

# build executable from Main.cpp
add_executable(test_main.exe Main.cpp)
target_link_libraries(test_main.exe PRIVATE fracLLL ${NTL_LIBRARY})
