project (solver_tools)

cmake_minimum_required (VERSION 3.14)

# Find eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_library (vector_tools REQUIRED)
find_library (error_tools REQUIRED)

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_subdirectory ("src/cpp")
add_subdirectory ("docs")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
