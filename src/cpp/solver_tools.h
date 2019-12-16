/*!
===============================================================================
|                               solver_tools.h                                |
===============================================================================
| The solver tools library. Incorporates a collection of non-linear solver    |
| tools built on top of Eigen. These tools are intended to be general enough  |
| to solve any small-ish nonlinear problem that can be solved using a Newton  |
| Raphson approach.                                                           |
===============================================================================
*/

#ifndef SOLVER_TOOLS_H
#define SOLVER_TOOLS_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>

namespace solverTools{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef double intType; //!Define the integer values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< intType > intVector; //!Define a vector of integers
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< intVector > intMatrix; //!Define a matrix of integers

    errorOut newtonRaphson( errorOut (*residual)(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                                                 floatVector &residual, floatMatrix &jacobian),
                            const floatVector &x0, 
                            floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, 
                            const unsigned int maxNLIterations=20, const floatType tolr=1e-9, const floatType tola=1e-9);

    errorOut checkTolerance( const floatVector &R, const floatVector &tol, bool &result);
}

#endif
