# solver\_tools

Tools for performing solves of nonlinear equations.

C20048 Tardigrade

 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

---

---

## Dependencies

### Executables

* CMake >= 3.14
* Doxygen >= 1.8.5

### Python Modules

For convenience, the minimal Python environment requirements for the
documentation and python interface builds are included in ``environment.yaml`` and
``requirements.txt``. A minimal anaconda environment for building the
documentation and python interface can be created from an existing anaconda
installation with the following commands.

```
$ conda env create --file environment.yaml
```

### Libraries

* eigen >= 3.3.7
* BOOST >= 1.53.0
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools

#### "Internal" project libraries

All of the ``{error,vector}_tools`` libraries are pulled from their git repos by
branch name and built with their respective cmake files as part of the cmake
build for this project.

#### Eigen

https://gitlab.com/libeigen/eigen

Eigen must be "installed" following the ``eigen/INSTALL`` instructions. The
Eigen dependence is easiest to resolve if eigen is installed in the default
install directory.  However, if you don't have admin privileges, you can also
insall Eigen to your home directory in ``~/include`` (or possibly in
``~/.local/include``, but this is untested by this project).

#### Non-admin Eigen install for solver_tools
[Reference](https://unix.stackexchange.com/questions/36871/where-should-a-local-executable-be-placed)

```
# Create personal include file directory
$ pwd
/home/$USER
$ mkdir .local/include
# Move to repository directory
$ cd /preferred/path/to/repos
# Example
$ pwd
# Clone eigen
$ git clone https://gitlab.com/libeigen/eigen.git
$ cd eigen
$ git checkout 3.3.7
# Build eigen
$ mkdir build
$ cd build
$ cmake3 .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
$ make install
```

---

---

## Build and Test

This repository is now built completely with cmake.

> **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
> API is much more useful

A build script has been created for convenience, ``new_build.sh``. It will build
everything including the library binary, the test binary, and the documentation.
This is the same build script used by ``build_test.sh`` for CI builds and
testing.

### build process

1) Activate the correct python environment

```
$ conda activate environment
```

2) Build everything

```
$ pwd
/path/to/solver_tools/

# Just perform the build
./new_build.sh

# Build and perform tests
./build_test.sh
```

3) View test results

```
cat build/src/cpp/tests/results.tex
```

4) Display docs

```
# Sphinx
firefox build/docs/sphinx/index.html &

# Doxygen
firefox build/docs/doxygen/html/index.html &
```

### Local development

In some cases it is not convenient to pull down every repository required but it may be desired that local
versions of the repository are used. An example of when this may be needed is if development is across
multiple libraries and is proceeding faster than collaborators can check in results. In this case, and
outside of developers no-one should need to do this, a version of the code using local repositories can be
built.

1) Activate a python environment

```
$ conda activate environment
```

2) Define convenience environment variables
$ my_error_tools=/path/to/my/error_tools
$ my_vector_tools=/path/to/my/vector_tools

3) Perform the initial configuration

```
$ pwd
/path/to/solver_tools
$ mkdir build
$ cd build
$ cmake .. -DFETCH_SOURCE=LOCAL -DERROR_TOOLS_PATH=${my_error_tools} -DVECTOR_TOOLS_PATH=${my_vector_tools}
```

4) Building the library

```
$ pwd
/path/to/solver_tools/build
$ make
```


### Building the documentation

To build just the documentation pick up the steps here:

2) Create the build directory and move there

```
$ pwd
/path/to/solver_tools/
$ mkdir build/
$ cd build/
```

3) Run cmake3 configuration

```
$ pwd
/path/to/solver_tools/build/
$ cmake3 ..
```

4) Build the docs

```
$ cmake3 --build docs
```

5) Documentation builds to:

```
solver_tools/build/docs/sphinx/index.html
```

6) Display docs

```
$ pwd
/path/to/solver_tools/build/
$ firefox docs/sphinx/index.html &
```

7) While the Sphinx API is still a WIP, try the doxygen API

```
$ pwd
/path/to/solver_tools/build/
$ firefox docs/doxygen/html/index.html &
```
