# - Try to find error_tools 
# Once done, this will define
#
#  ERROR_TOOLS_FOUND - system has ERROR_TOOLS
#  ERROR_TOOLS_INCLUDE_DIRS - the ERROR_TOOLS include directories
#  ERROR_TOOLS_LIBRARIES - link these to use ERROR_TOOLS

include(LibFindMacros)

# Dependencies
libfind_package(error_tools)

# Set hard coded relative path
set(ERROR_TOOLS_RELATIVE_PATH "${PROJECT_SOURCE_DIR}/../error_tools/src/cpp")

# Include dir
find_path(ERROR_TOOLS_INCLUDE_DIR
          NAMES error_tools.h
          PATHS ${ERROR_TOOLS_RELATIVE_PATH})

find_library(ERROR_TOOLS_LIBRARY
             NAMES error_tools
             PATHS ${ERROR_TOOLS_RELATIVE_PATH})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ERROR_TOOLS_PROCESS_INCLUDES ERROR_TOOLS_INCLUDE_DIR)
set(ERROR_TOOLS_PROCESS_LIBS ERROR_TOOLS_LIBRARY)
libfind_process(ERROR_TOOLS)
