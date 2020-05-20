# - Try to find vector_tools 
# Once done, this will define
#
#  VECTOR_TOOLS_FOUND - system has VECTOR_TOOLS
#  VECTOR_TOOLS_INCLUDE_DIRS - the VECTOR_TOOLS include directories
#  VECTOR_TOOLS_LIBRARIES - link these to use VECTOR_TOOLS

include(LibFindMacros)

# Dependencies
libfind_package(vector_tools error_tools)

# Set hard coded relative path
set(VECTOR_TOOLS_RELATIVE_PATH "${PROJECT_SOURCE_DIR}/../vector_tools/src/cpp")

# Include dir
find_path(VECTOR_TOOLS_INCLUDE_DIR
          NAMES vector_tools.h
          PATHS ${VECTOR_TOOLS_RELATIVE_PATH})

find_library(VECTOR_TOOLS_LIBRARY
             NAMES vector_tools
             PATHS ${VECTOR_TOOLS_RELATIVE_PATH})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(VECTOR_TOOLS_PROCESS_INCLUDES VECTOR_TOOLS_INCLUDE_DIR ERROR_TOOLS_INCLUDE_DIRS)
set(VECTOR_TOOLS_PROCESS_LIBS VECTOR_TOOLS_LIBRARY ERROR_TOOLS_LIBRARIES)
libfind_process(VECTOR_TOOLS)
