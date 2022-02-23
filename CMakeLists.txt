cmake_minimum_required(VERSION 3.11...3.22)

project(
    fortune-hyperbolic-sample-application
    VERSION 1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

add_subdirectory(generator)
add_subdirectory(experiments)

add_executable(main main.cpp)

set(FORTUNE_USE_MPFR TRUE)
list(APPEND CMAKE_PREFIX_PATH fortune-hyperbolic)
find_package(Fortune_Hyperbolic)

target_include_directories(main PRIVATE ${Fortune_Hyperbolic_INCLUDE_DIRS})
target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/include)

message(STATUS ${Fortune_Hyperbolic_INCLUDE_DIRS})
message(STATUS ${Fortune_Hyperbolic_LIBRARIES})

target_link_libraries(main ${Fortune_Hyperbolic_LIBRARIES})
target_compile_options(main PUBLIC -Wall)

add_subdirectory(tests)