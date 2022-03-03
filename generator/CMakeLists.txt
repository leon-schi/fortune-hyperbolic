add_executable(generator generator.cpp)
target_include_directories(generator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_compile_options(generator PUBLIC -Wall)