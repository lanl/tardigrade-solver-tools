project (solver_tools)

cmake_minimum_required (VERSION 3.14)

# Find eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Find related, but external, projects
include(FetchContent)
FetchContent_Declare(
    error_tools
#    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/error_tools.git
#    GIT_TAG        origin/dev
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/~kbrindley/error_tools.git
    GIT_TAG        origin/feature/VIP-293-error_tools-make-docs-and-tests-optional-build-options
)
FetchContent_Declare(
    vector_tools
#    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/vector_tools.git
#    GIT_TAG        origin/dev
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/~kbrindley/vector_tools.git
    GIT_TAG        origin/feature/VIP-294-vector_tools-make-docs-and-tests-optional-build-options
)
FetchContent_MakeAvailable(error_tools vector_tools)

# Add project source directories
add_subdirectory ("src/cpp")
add_subdirectory ("docs")

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
