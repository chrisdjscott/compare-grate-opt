cmake_minimum_required(VERSION 3.1.0)
project(testopt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# default to RELEASE build if not specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

find_package(GSL REQUIRED)

add_executable(testopt test.cpp)
target_link_libraries(testopt GSL::gsl)
