FILE(GLOB HEADER_LIST "${fortune-hyperbolic_SOURCE_DIR}/include/*.hpp")
FILE(GLOB SOURCE_LIST "${fortune-hyperbolic_SOURCE_DIR}/src/*.cpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(fortune_hyperbolic ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(fortune_hyperbolic PUBLIC "${PROJECT_BINARY_DIR}/include")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
find_package(Boost 1.73 REQUIRED)
#find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
#include(../FindMPFR.cmake)
target_include_directories(fortune_hyperbolic PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(fortune_hyperbolic LINK_PUBLIC ${MPFR_LIBRARIES})

target_compile_features(fortune_hyperbolic PUBLIC cxx_std_17)
target_compile_options(fortune_hyperbolic PRIVATE -Wall)
