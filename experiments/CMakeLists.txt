set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
include(CGAL-5.2/CGALConfig.cmake)

add_executable(old old.cpp)
target_include_directories(old PRIVATE ./CGAL-5.2/include)
target_link_libraries(old CGAL)