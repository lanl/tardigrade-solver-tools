add_library (solver_tools "solver_tools.cpp" "solver_tools.h")

target_include_directories(solver_tools PUBLIC .)
set_target_properties(solver_tools PROPERTIES CXX_STANDARD 11)
target_compile_options(solver_tools PUBLIC
                       -Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb)

# FIXME: Hack in the include directories until find_packages works
set(VECTOR_TOOLS_ABS_PATH
    /Users/projects/kbrindley/e13repos/constitutive_models/vector_tools/src/cpp)
set(ERROR_TOOLS_ABS_PATH
    /Users/projects/kbrindley/e13repos/constitutive_models/error_tools/src/cpp)

target_link_libraries(solver_tools Eigen3::Eigen
                      "${ERROR_TOOLS_ABS_PATH}/error_tools.o")
# TODO: Replace include directories with fing_package paths when they work
#                      "${VECTOR_TOOLS_LIBRARIES}" "${ERROR_TOOLS_LIBRARIES}")

target_include_directories(solver_tools PUBLIC "${VECTOR_TOOLS_ABS_PATH}" "${ERROR_TOOLS_ABS_PATH}")
# TODO: Replace include directories with fing_package paths when they work
#                           "${VECTOR_TOOLS_INCLUDE_DIRS}"
#                           "${ERROR_TOOLS_INCLUDE_DIRS}")
