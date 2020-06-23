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
#    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/~kbrindley/error_tools.git
#    GIT_TAG        origin/feature/VIP-293-error_tools-make-docs-and-tests-optional-build-options
    GIT_REPOSITORY /Users/projects/kbrindley/e13repos/constitutive_models/error_tools
    GIT_TAG        feature/VIP-293-error_tools-make-docs-and-tests-optional-build-options
)
FetchContent_Declare(
    vector_tools
#    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/vector_tools.git
#    GIT_TAG        origin/dev
#    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/~kbrindley/vector_tools.git
#    GIT_TAG        origin/feature/VIP-294-vector_tools-make-docs-and-tests-optional-build-options
    GIT_REPOSITORY /Users/projects/kbrindley/e13repos/constitutive_models/vector_tools
    GIT_TAG        feature/VIP-294-vector_tools-make-docs-and-tests-optional-build-options
)
FetchContent_MakeAvailable(error_tools vector_tools)
FetchContent_GetProperties(error_tools)
if(NOT error_tools_POPULATED)
    FetchContent_Populate(error_tools)
    add_subdirectory(${error_tools_SOURCE_DIR})
endif()
FetchContent_GetProperties(vector_tools)
if(NOT vector_tools_POPULATED)
    FetchContent_Populate(vector_tools)
    add_subdirectory(${vector_tools_SOURCE_DIR})
endif()

# Add project source directories
add_subdirectory ("src/cpp")
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory ("docs")
endif()

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
