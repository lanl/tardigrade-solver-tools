# solver\_tools

Tools for performing solves of nonlinear equations. 

Note: In order to use the Intel compiler one must run the following command 
in a bash prompt:
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

This is the same command that the abaqus command issues. It may be that 
this command will change on different platforms.

---

---

## Dependencies

### Make

These tools have several dependencies that must be available in the same parent
directory as this repo. 

* eigen: https://gitlab.com/libeigen/eigen
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools

### CMake

The project is transitioning from Make to Cmake. For cmake builds, Eigen must be
"installed" following the ``eigen/INSTALL`` instructions. The Eigen dependence
is easiest to resolve if eigen is installed in the default install directory.

---

---

## Documentation

The documentation for this project is built with cmake, sphinx, doxygen, and
breathe.

The documentation was built with these [microsoft developer blog
instructions](https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/)
for c++ projects. The one caveat is the ``vector_tools`` project is used as a
header only project, so no libraries are built.

To build the documentation run the following from the project root directory,
``vector_tools``:

```
$ pwd
path/to/vector_tools
$ mkdir build
$ cd build
$ cmake3 ..
$ cmake --build docs
```
