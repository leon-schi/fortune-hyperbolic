set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
include(CGAL-5.2/CGALConfig.cmake)

add_executable(old old.cpp)
target_include_directories(old PRIVATE ${CGAL_INCLUDE_DIRS})
target_include_directories(old PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(old CGAL)