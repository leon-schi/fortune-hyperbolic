FILE(GLOB HEADER_LIST "${fortune-hyperbolic_SOURCE_DIR}/include/*.hpp")
FILE(GLOB SOURCE_LIST "${fortune-hyperbolic_SOURCE_DIR}/src/*.cpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(fortune_hyperbolic ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(fortune_hyperbolic PUBLIC ../include)
target_include_directories(fortune_hyperbolic PUBLIC "${PROJECT_BINARY_DIR}/include")

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(fortune_hyperbolic PUBLIC
                -Wall -Wextra -Wfloat-equal -Wconversion
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
                -fno-sanitize-recover=all
        )
        target_link_options(fortune_hyperbolic PUBLIC
                -fsanitize=address,undefined
                -fno-omit-frame-pointer)
    endif()
endif()

target_compile_features(fortune_hyperbolic PUBLIC cxx_std_17)
target_compile_options(fortune_hyperbolic PRIVATE -Wall)
