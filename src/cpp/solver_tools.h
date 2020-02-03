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
#include<sstream>

namespace solverTools{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef int intType; //!Define the integer values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< intType > intVector; //!Define a vector of integers
    typedef std::vector< floatVector > floatMatrix; //!Define a matrix of floats
    typedef std::vector< intVector > intMatrix; //!Define a matrix of integers

    typedef errorOut(*NonLinearFunction)(const floatVector&, const floatMatrix&, const intMatrix&, floatVector&);
    typedef std::function< errorOut(const floatVector&, const floatMatrix&, const intMatrix&, floatVector&) > stdFncNLF;
    typedef errorOut(*NonLinearFunctionWithJacobian)(const floatVector&, const floatMatrix&, const intMatrix&, floatVector&, floatMatrix&, 
                                                     floatMatrix &, intMatrix &);
    typedef std::function< errorOut(const floatVector&, const floatMatrix&, const intMatrix&, floatVector&, 
                                    floatMatrix&, floatMatrix&, intMatrix&) > stdFncNLFJ;

    errorOut newtonRaphson( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                            const floatVector &x0, 
                            floatVector &x, bool &convergeFlag, floatMatrix &floatOuts, intMatrix &intOuts, 
                            const floatMatrix &floatArgs, const intMatrix &intArgs, 
                            const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                            const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5);

    errorOut homotopySolver( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                            const floatVector &x0,
                            floatVector &x, bool &convergeFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                            const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const unsigned int homotopySteps=10);

    errorOut checkTolerance( const floatVector &R, const floatVector &tol, bool &result);

    errorOut checkLSCriteria( const floatVector &R, const floatVector &Rp, bool &result, const floatType alpha=1e-4);

    errorOut finiteDifference( stdFncNLF fxn,
                            const floatVector &x0,
                            floatMatrix &J, const floatMatrix &floatArgs, const intMatrix &intArgs, const floatType eps=1e-6);

    errorOut checkJacobian( stdFncNLFJ residual,
                            const floatVector &x0,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, bool &isGood, const floatType eps=1e-6,
                            const floatType tolr=1e-6, const floatType tola=1e-6, const bool suppressOutput = false);

}

#endif
