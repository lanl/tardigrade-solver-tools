/**
  *****************************************************************************
  * \file solver_tools.h
  *****************************************************************************
  * The solver tools library. Incorporates a collection of non-linear solver
  * tools built on top of Eigen. These tools are intended to be general enough
  * to solve any small-ish nonlinear problem that can be solved using a Newton
  * Raphson approach.
  *****************************************************************************
  */

#ifndef SOLVER_TOOLS_H
#define SOLVER_TOOLS_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<sstream>

namespace solverTools{

    typedef errorTools::Node errorNode; //!< Redefinition for the error node
    typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
    typedef double floatType; //!< Define the float values type.
    typedef int intType; //!< Define the integer values type.
    typedef std::vector< floatType > floatVector; //!<  Define a vector of floats
    typedef std::vector< intType > intVector; //!< Define a vector of integers
    typedef std::vector< floatVector > floatMatrix; //!< Define a matrix of floats
    typedef std::vector< intVector > intMatrix; //!< Define a matrix of integers

    using solverType = vectorTools::solverType< floatType >; //!< Force consistency with vectorTools

    /**
     * A residual function.
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &residual: The residual vector
     */
    typedef errorOut(*NonLinearFunction)(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, floatVector &residual);
    /**
     * A residual function
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &residual: The residual vector
     */
    typedef std::function< errorOut(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, floatVector &residual) > stdFncNLF;

    /**
     * A residual function including the Jacobian
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &residual: The residual vector
     * \param &jacobian: The jacobian matrix of the residual w.r.t. x
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */
    typedef errorOut(*NonLinearFunctionWithJacobian)(const floatVector&, const floatMatrix&, const intMatrix&, floatVector&, floatMatrix&,
                                                     floatMatrix &, intMatrix &);
    /**
     * A residual function including the Jacobian
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &residual: The residual vector
     * \param &jacobian: The jacobian matrix of the residual w.r.t. x
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */
    typedef std::function< errorOut(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, floatVector &residual,
                                    floatMatrix &jacobian, floatMatrix &floatOuts, intMatrix &intOuts) > stdFncNLFJ;
    /**
     * A Lagrangian which also returns its gradient.
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param intMatrix &intArgs: Additional integer arguments to the residual
     * \param &value: The value of the Lagrangian.
     * \param &gradient: The gradient of the Lagrangian function.
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */
    typedef errorOut(*LagrangianFunctionWithGradient)(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, floatType &value, floatVector &gradient,
                                                     floatMatrix &floatOuts, intMatrix &intOuts);

    /**
     * A Lagrangian which also returns its gradient.
     * 
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param intMatrix &intArgs: Additional integer arguments to the residual
     * \param &value: The value of the Lagrangian.
     * \param &gradient: The gradient of the Lagrangian function.
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */
    typedef std::function< errorOut(const floatVector&, const floatMatrix&, const intMatrix&, floatType&, floatVector&,
                                    floatMatrix&, intMatrix&) > stdFncLagrangianG;

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                            const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false );

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                            const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                            const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false );

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                            const intVector &boundVariableIndices, const intVector &boundSigns, const floatVector &boundValues,
                            const bool boundMode,
                            const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                            const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false );

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs,
                             const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                             const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const floatType ds0 = 1.0,
                             const floatType dsMin = 0.1, const bool resetOuts = false );

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                             const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                             const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const floatType ds0 = 1.0,
                             const floatType dsMin = 0.1, const bool resetOuts = false );

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                             const intVector &boundVariableIndices, const intVector &boundSigns, const floatVector &boundValues,
                             const bool boundMode,
                             const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                             const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const floatType ds0 = 1.0,
                             const floatType dsMin = 0.1, const bool resetOuts = false );

    errorOut barrierHomotopySolver( stdFncNLFJ residual, const floatType &dt, const floatVector &x0,
                                    const intVector &variableIndices, const intVector &residualIndices,
                                    const intVector &barrierSigns, const floatVector &barrierValues,
                                    const floatVector &logAMaxValues,
                                    const floatMatrix &floatArgs, const intMatrix &intArgs,
                                    const bool &implicitRefine,
                                    floatVector &x, bool &convergeFlag, bool &fatalErrorFlag,
                                    floatMatrix &floatOuts, intMatrix &intOuts,
                                    const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                                    const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false );

    errorOut barrierHomotopySolver( stdFncNLFJ residual, const floatType &dt, const floatVector &x0,
                                    const intVector &variableIndices, const intVector &residualIndices,
                                    const intVector &barrierSigns, const floatVector &barrierValues,
                                    const floatVector &logAMaxValues,
                                    const floatMatrix &floatArgs, const intMatrix &intArgs,
                                    const bool &implicitRefine,
                                    floatVector &x, bool &convergeFlag, bool &fatalErrorFlag,
                                    floatMatrix &floatOuts, intMatrix &intOuts, solverType &linearSolver, floatMatrix &J,
                                    const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                                    const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false );

    errorOut BFGS( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                           floatType &, floatVector &, floatMatrix &, intMatrix &
                                         ) > lagrangianGradient,
                   const floatVector &x0,
                   floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                   const floatMatrix &floatArgs, const intMatrix &intArgs,
                   const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                   const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const bool resetOuts = false,
                   const floatType stepSize = 1., const floatType maxdx = -1 );

    errorOut homotopyBFGS( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatType &, floatVector &, floatMatrix &, intMatrix &
                                                   ) > lagrangianGradientFunction,
                           const floatVector &x0,
                           floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                           const floatMatrix &floatArgs, const intMatrix &intArgs,
                           const unsigned int maxNLIterations = 20, const floatType tolr = 1e-9, const floatType tola = 1e-9,
                           const floatType alpha = 1e-4, const unsigned int maxLSIterations = 5, const floatType ds0 = 1.,
                           const floatType dsMin = 0.1, const bool resetOuts = false, const floatType maxdx = -1 );

    errorOut checkTolerance( const floatVector &R, const floatVector &tol, bool &result);

    errorOut checkLSCriteria( const floatVector &R, const floatVector &Rp, bool &result, const floatType alpha=1e-4);

    errorOut finiteDifference( stdFncNLF fxn,
                            const floatVector &x0,
                            floatMatrix &J, const floatMatrix &floatArgs, const intMatrix &intArgs, const floatType eps=1e-6);

    errorOut checkJacobian( stdFncNLFJ residual,
                            const floatVector &x0,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, bool &isGood, const floatType eps=1e-6,
                            const floatType tolr=1e-6, const floatType tola=1e-6, const bool suppressOutput = false);

    errorOut aFxn( const floatType &pseudoT, const floatType logAmax, floatType &a );

    errorOut aFxn( const floatType &pseudoT, const floatType logAMax, floatType &a, floatType &dadt );

    errorOut computeBarrierFunction( const floatType &x, const floatType &pseudoT, const floatType &logAmax,
                                     const floatType &b, const bool &sign, floatType &barrierFunction );

    errorOut computeBarrierFunction( const floatType &x, const floatType &pseudoT, const floatType &logAmax,
                                     const floatType &b, const bool &sign, floatType &barrierFunction,
                                     floatType &dbdx, floatType &dbdt );

    errorOut computeBarrierHomotopyResidual( stdFncNLFJ computeOriginalResidual, const floatVector &x,
                                             const floatMatrix &floatArgs, const intMatrix &intArgs,
                                             floatVector &residual, floatMatrix &jacobian,
                                             floatMatrix &floatOuts, intMatrix &intOuts
                                           );

    errorOut applyBoundaryLimitation( const floatVector &x0, const intVector &variableIndices, const intVector &barrierSigns,
                                      const floatVector &barrierValues, floatVector &dx,
                                      const floatType tolr = 1e-9, const floatType tola = 1e-9, const bool mode = false );

}

#endif
