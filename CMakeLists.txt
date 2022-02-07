cmake_minimum_required(VERSION 3.11...3.22)

project(
    fortune-hyperbolic
    VERSION 1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(include/cmake.hpp.in include/cmake.hpp)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
                -Wall -Wextra -Wfloat-equal -Wconversion
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
                -fno-sanitize-recover=all)
        add_link_options(
                -fsanitize=address,undefined
                -fno-omit-frame-pointer)
    endif()
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(generator)
add_subdirectory(experiments)

add_executable(main main.cpp)

target_link_libraries(main fortune_hyperbolic)
target_compile_options(main PUBLIC -Wall)