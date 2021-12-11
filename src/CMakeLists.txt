set(HEADER_LIST "${fortune-hyperbolic_SOURCE_DIR}/include/lib.hpp")

add_library(fortune_hyperbolic lib.cpp ${HEADER_LIST})
target_include_directories(fortune_hyperbolic PUBLIC ../include)
target_include_directories(fortune_hyperbolic PUBLIC "${PROJECT_BINARY_DIR}/include")

target_compile_features(fortune_hyperbolic PUBLIC cxx_std_17)
target_compile_options(fortune_hyperbolic PRIVATE -Wall)
