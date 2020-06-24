add_library (solver_tools "solver_tools.cpp" "solver_tools.h")

target_include_directories(solver_tools PUBLIC .)
set_target_properties(solver_tools PROPERTIES CXX_STANDARD 11)
target_compile_options(solver_tools PUBLIC
                       -Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb)

target_link_libraries(solver_tools Eigen3::Eigen)
target_include_directories(solver_tools PUBLIC "${error_tools_SOURCE_DIR}/src/cpp" "${vector_tools_SOURCE_DIR}/src/cpp")
