cmake_minimum_required(VERSION 3.11...3.22)

project(
        fortune-hyperbolic
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

get_filename_component(FORTUNE_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

configure_file(${FORTUNE_CONFIG_DIR}/include/fortune-hyperbolic/cmake.hpp.in include/cmake.hpp)

list(APPEND Fortune_Hyperbolic_INCLUDE_DIRS ${FORTUNE_CONFIG_DIR}/include)
list(APPEND Fortune_Hyperbolic_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include)

if (FORTUNE_USE_MPFR)
    find_package(Boost 1.73 REQUIRED)
    include(${FORTUNE_CONFIG_DIR}/FindMPFR.cmake)
    list(APPEND Fortune_Hyperbolic_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
    set(Fortune_Hyperbolic_LIBRARIES ${MPFR_LIBRARIES})
endif()