add_library (solver_tools "solver_tools.cpp" "solver_tools.h")

target_include_directories(solver_tools PUBLIC .)
set_target_properties(solver_tools PROPERTIES CXX_STANDARD 11)

target_link_libraries (solver_tools Eigen3::Eigen)

target_compile_options(solver_tools PUBLIC
                       -Wall -ansi -pedantic -O3 -fmax-errors=5)
