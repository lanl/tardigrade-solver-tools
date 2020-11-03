/**
  *****************************************************************************
  * \file solver_tools.cpp
  *****************************************************************************
  * The solver tools library. Incorporates a collection of non-linear solver
  * tools built on top of Eigen. These tools are intended to be general enough
  * to solve any small-ish nonlinear problem that can be solved using a Newton
  * Raphson approach.
  *****************************************************************************
  */

#include<solver_tools.h>

namespace solverTools{

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation
         * of a typical Newton-Raphson solver which can take an arbitrary
         * residual function.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if a fatal error has occurred
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param resetOuts: A flag for whether the output matrices should be reset
         *     prior to each iteration.
         */

        solverType linearSolver;
        floatMatrix J;
        return newtonRaphson( residual, x0, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                              maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts );
    }

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation
         * of a typical Newton-Raphson solver which can take an arbitrary
         * residual function.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if a fatal error has occurred
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param &linearSolver: The linear solver used in the solve. This object
         *     contains the decomposed Jacobian matrix which can be very useful in 
         *     total Jacobian calculations
         * \param &J: The Jacobian matrix. Useful in calculating total Jacobians.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param resetOuts: A flag for whether the output matrices should be reset
         *     prior to each iteration.
         */

        intVector boundVariableIndices(0);
        intVector boundSigns(0);
        floatVector boundValues(0);
        bool boundMode = false;
        return newtonRaphson( residual, x0, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                              boundVariableIndices, boundSigns, boundValues, boundMode,
                              maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts );
    }

    errorOut newtonRaphson( stdFncNLFJ residual, const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                            const intVector &boundVariableIndices, const intVector &boundSigns, const floatVector &boundValues,
                            const bool boundMode,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation
         * of a typical Newton-Raphson solver which can take an arbitrary
         * residual function.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if a fatal error has occurred
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param &linearSolver: The linear solver used in the solve. This object
         *     contains the decomposed Jacobian matrix which can be very useful in 
         *     total Jacobian calculations
         * \param &J: The Jacobian matrix. Useful in calculating total Jacobians
         * \param &boundVariableIndices: The indices of variables that have hard bounds
         * \param &boundSigns: The signs of the bounds. 0 for a negative ( lower ) bound
         *     and 1 for a positive ( upper ) bound
         * \param &boundValues: The values of the boundaries.
         * \param boundMode: The mode for the boundary. See applyBoundaryLimitation for
         *     more details.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param resetOuts: A flag for whether the output matrices should be reset
         *     prior to each iteration
         */

        //Compute the initial residual and jacobian
        floatVector dx = floatVector(x0.size(), 0);
        floatVector ddx;
        floatVector R, Rp;

        //Make copies of the initial float out values
        floatMatrix oldFloatOuts = floatOuts;
        intMatrix   oldIntOuts   = intOuts;

        errorOut error = residual(x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts);

        if ( error ){
            if ( ( R.size( ) == 101 ) && ( J.size( ) == 212 ) ){ // This is supposed to detect a failure in convergence of a sub non-linear solve. I don't like this approach and will be setting up an issue to change it
                convergeFlag = false;
                fatalErrorFlag = false;
                errorOut result = new errorNode( "newtonRaphson", "Convergence error in intial residual" );
                result->addNext( error );
                return result;
            }
            else {
                fatalErrorFlag = true;
                errorOut result = new errorNode("newtonRaphson", "Error in computation of initial residual");
                result->addNext(error);
                return result;
            }
        }

        if (R.size() != x0.size()){
            fatalErrorFlag = true;
            return new errorNode("newtonRaphson", "The residual and x0 don't have the same lengths. The problem is ill-defined.");
        }

        if ( resetOuts ){
            floatOuts = oldFloatOuts;
            intOuts = oldIntOuts;
        }

        //Set the tolerance for each value individually
        floatVector tol = floatVector( R.size( ), 0 );
        for ( unsigned int i = 0; i < R.size( ); i++ ){ tol[ i ] = tolr * fabs( R[ i ] ) + tola; }

        //Copy R to Rp
        Rp = R;

        //Initialize variables required for the iteration loop
        unsigned int nNLIterations = 0;
        unsigned int nLSIterations = 0;
        float lambda = 1;
        bool converged, lsCheck;
        checkTolerance( R, tol, converged );
        unsigned int rank;
        convergeFlag = false;
        fatalErrorFlag = false;

        //Begin the iteration loop
        while ( ( !converged ) && ( nNLIterations<maxNLIterations ) ){

            //Perform the linear solve
            ddx = -vectorTools::solveLinearSystem( J, R, rank, linearSolver );

            //Check the rank to make sure the linear system has a unique solution
            if ( rank != R.size( ) ){
                convergeFlag = false;
                return new errorNode( "newtonRaphson", "The jacobian matrix is singular" );
            }

            //Apply any boundaries on the variables
            error = applyBoundaryLimitation( x0 + dx, boundVariableIndices, boundSigns, boundValues, ddx,
                                             tolr, tola, boundMode );

            if ( error ){
                errorOut result = new errorNode( "newtonRaphson", "Error in the application of the boundary limitations" );
                result->addNext( error );
                return result;
            }

            //Update dx
            dx += ddx;

            if ( resetOuts ){
                floatOuts = oldFloatOuts;
                intOuts = oldIntOuts;
            }

            //Compute the new residual
            error = residual(x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts);

            if (error){
                if ( ( R.size( ) == 101 ) && ( J.size( ) == 212 ) ){ //TODO: Replace with something better
                    errorOut result = new errorNode( "newtonRaphson", "Convergence error in sub Newton-Raphson process" );
                    result->addNext( error );
                    fatalErrorFlag = false;
                    convergeFlag = false;
                    return result;
                }
                else{
                    fatalErrorFlag = true;
                    errorOut result = new errorNode("newtonRaphson", "Error in residual calculation in non-linear iteration");
                    result->addNext(error);
                    return result;
                }
            }

            //Check the line search criteria
            checkLSCriteria( R, Rp, lsCheck, alpha );
            nLSIterations = 0;
            lambda = 1;

            //Enter line-search if required
            while ( ( !lsCheck ) && ( nLSIterations < maxLSIterations ) ){

                //Extract ddx from dx
                dx -= lambda * ddx;

                //Decrement lambda. We could make this fancier but just halving it is probably okay
                lambda *= 0.5;

                //Update dx
                dx += lambda * ddx;

                //Reset floatOuts and intOuts to the previously converged values
                floatOuts = oldFloatOuts;
                intOuts   = oldIntOuts;

                //Compute the new residual
                error = residual( x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts );

                if ( error ){
                    if ( ( R.size( ) == 101 ) && ( J.size( ) == 212 ) ){//TODO: I continue to hate this
                        errorOut result = new errorNode( "newtonRaphson", "Convergence error in residual function" );
                        result->addNext( error );
                        fatalErrorFlag = false;
                        convergeFlag = false;
                        return result;
                    }
                    else{
                        fatalErrorFlag = true;
                        errorOut result = new errorNode( "newtonRaphson", "Error in line-search" );
                        result->addNext( error );
                        return result;
                    }
                }

                //Perform the line-search check
                checkLSCriteria( R, Rp, lsCheck, alpha );

                //Increment the number of line-search iterations
                nLSIterations++;
            }

            if ( !lsCheck ){
                convergeFlag = false;
                fatalErrorFlag = false;
                return new errorNode( "newtonRaphson", "The line-search failed to converge." );
            }
            else{
                Rp = R;
                if ( !resetOuts ){
                    oldFloatOuts = floatOuts;
                    oldIntOuts   = intOuts;
                }
            }

            //Check if the solution is converged
            checkTolerance( R, tol, converged );

            //Increment nNLIterations
            nNLIterations++;
        }

        //Check if the solution converged
        if ( !converged ){
            convergeFlag = false;
            return new errorNode( "newtonRaphson", "The Newton-Raphson solver failed to converge." );
        }
        else{
            //Update x
            x = x0 + dx;
            //Solver completed successfully
            convergeFlag = true;
            return NULL;
        }
    }

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs,
                             const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                             const floatType alpha, const unsigned int maxLSIterations,const floatType ds0,
                             const floatType dsMin, const bool resetOuts ){
        /*!
         * Solve a non-linear equation using a homotopy Newton solver. This method
         * can be successful in solving very stiff equations which other techniques
         * struggle to capture. It effectively breaks the solve up into sub-steps
         * of easier to solve equations which will eventually converge to the
         * more difficult problem.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if there has been a fatal error in the solve
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param ds0: The initial step size.
         * \param dsMin: The minimum step size.
         * \param resetOuts: Flag for whether the outputs should be reset at each step
         */

        solverType linearSolver;
        floatMatrix J;
        return homotopySolver( residual, x0, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                               maxNLIterations, tolr, tola, alpha, maxLSIterations, ds0, dsMin, resetOuts );

    }

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                             const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                             const floatType alpha, const unsigned int maxLSIterations,const floatType ds0,
                             const floatType dsMin, const bool resetOuts ){
        /*!
         * Solve a non-linear equation using a homotopy Newton solver. This method
         * can be successful in solving very stiff equations which other techniques
         * struggle to capture. It effectively breaks the solve up into sub-steps
         * of easier to solve equations which will eventually converge to the
         * more difficult problem.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function.
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if there has been a fatal error in the solve
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param &linearSolver: The linear solver object used in the solution of the nonlinear solve.
         *     Note that the linear solver will always be the same as that is defined in solverTools.
         *     This contains the decomposed matrix which is useful for Jacobian calculations.
         * \param &J: The Jacobian matrix of the nonlinear solve.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param ds0: The initial step size.
         * \param dsMin: The minimum step size.
         * \param resetOuts: Flag for whether the outputs should be reset at each step
         */

        intVector boundVariableIndices(0);
        intVector boundSigns(0);
        floatVector boundValues(0);
        bool boundMode = false;
        return homotopySolver( residual, x0, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts, floatArgs, intArgs, linearSolver, J,
                               boundVariableIndices, boundSigns, boundValues, boundMode,
                               maxNLIterations, tolr, tola, alpha, maxLSIterations, ds0, dsMin, resetOuts );

    }

    errorOut homotopySolver( stdFncNLFJ residual, const floatVector &x0,
                             floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                             const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                             const intVector &boundVariableIndices, const intVector &boundSigns, const floatVector &boundValues,
                             const bool boundMode,
                             const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                             const floatType alpha, const unsigned int maxLSIterations, const floatType ds0,
                             const floatType dsMin, const bool resetOuts ){
        /*!
         * Solve a non-linear equation using a homotopy Newton solver. This method
         * can be successful in solving very stiff equations which other techniques
         * struggle to capture. It effectively breaks the solve up into sub-steps
         * of easier to solve equations which will eventually converge to the
         * more difficult problem.
         *
         * The main routine accepts the following parameters:
         * \param residual: The residual function.
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &fatalErrorFlag: A flag which indicates if there has been a fatal error in the solve
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param &linearSolver: The linear solver object used in the solution of the nonlinear solve.
         *     Note that the linear solver will always be the same as that is defined in solverTools.
         *     This contains the decomposed matrix which is useful for Jacobian calculations.
         * \param &J: The Jacobian matrix of the nonlinear solve.
         * \param &boundVariableIndices: The indices of variables that have hard bounds
         * \param &boundSigns: The signs of the bounds. 0 for a negative ( lower ) bound
         *     and 1 for a positive ( upper ) bound
         * \param &boundValues: The values of the boundaries.
         * \param boundMode: The mode for the boundary. See applyBoundaryLimitation for
         *     more details.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param ds0: The initial step size.
         * \param dsMin: The minimum step size.
         * \param resetOuts: Flag for whether the outputs should be reset at each step
         */

        //Initialize the homotopy solver
        floatType ds = ds0;
        floatType s  = 0;
        floatVector xh = x0;
        floatVector xdot = x0;
        floatVector rh, dxh;
        solverType ls;
        unsigned int rank;

        //Save the floatOuts and intOuts
        floatMatrix oldFloatOuts = floatOuts;
        intMatrix   oldIntOuts   = intOuts;

        //Define the homotopy residual equation
        floatVector Rinit;

        //Compute the initial residual
        errorOut error = residual( x0, floatArgs, intArgs, Rinit, J, floatOuts, intOuts );

        if ( error ){
            fatalErrorFlag = true;
            errorOut result = new errorNode( "homotopySolver", "error in initial residual calculation" );
            result->addNext( error );
            return result;
        }

        //Define the homotopy residual
        stdFncNLFJ homotopyResidual;
        homotopyResidual = [&]( const floatVector &homotopy_x, const floatMatrix &homotopy_floatArgs, const intMatrix &homotopy_intArgs,
                                floatVector &homotopy_residual, floatMatrix &homotopy_J, floatMatrix &homotopy_floatOuts, intMatrix &homotopy_intOuts ){
            /*!
             * A sub-function that computes the homotopy residual. This residual takes the original function and maps
             * it to another that is easier to solve. The variables into this function are all the same as in the
             * original residual but we modify the residual to include an offset which makes the expression
             * easier to solve.
             */

            floatVector R;

            error = residual( homotopy_x, homotopy_floatArgs, homotopy_intArgs, R, homotopy_J, homotopy_floatOuts, homotopy_intOuts );

            if ( error ){
                homotopy_residual = R;
                errorOut result = new errorNode( "homotopySolver::homotopyResidual", "error in residual calculation" );
                result->addNext( error );
                return result;
            }

            homotopy_residual = R - ( 1 - s ) * Rinit;

            return static_cast< errorOut >( NULL );
        };

        //Begin the homotopy loop
        while ( s < 1 ){

            //Update s
            s += ds;
            s = std::min( s, 1. );

            //Initialize the solver
            convergeFlag = false;

            if ( !resetOuts ){
                oldFloatOuts = floatOuts;
                oldIntOuts   = intOuts;
            }
            else{
                floatOuts = oldFloatOuts;
                intOuts   = oldIntOuts;
            }

            //Begin the adaptive homotopy loop
            while ( !convergeFlag ){

                //Compute the explicit estimate of xh ( this is kind of what makes it a homotopy )
                error = homotopyResidual( xh, floatArgs, intArgs, rh, J, floatOuts, intOuts );
    
                if ( error ){
                    fatalErrorFlag = true;
                    convergeFlag = false;
                    errorOut result = new errorNode( "homotopySolver", "The explicit homotopy estimate of x failed in an unexpected way. This shouldn't happen." );
                    result->addNext( error );
                    return result;
                }

                floatOuts = oldFloatOuts;
                intOuts   = oldIntOuts;
    
                xdot = -vectorTools::solveLinearSystem( J, rh, rank );
   
                if ( rank != J.size( ) ){
                    convergeFlag = false;
                    fatalErrorFlag = false;
                }
                else{
    
                    dxh = ds * xdot;

                    error = applyBoundaryLimitation( xh, boundVariableIndices, boundSigns, boundValues, dxh, tolr, tola, boundMode );
    
                    if ( error ){
                        errorOut result = new errorNode( "homotopySolve", "Fatal error in application of boundary limitations" );
                        result->addNext( error );
                        fatalErrorFlag = true;
                        return result;
                    }

                    error = newtonRaphson( homotopyResidual, xh + dxh, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts,
                                           floatArgs, intArgs, linearSolver, J, boundVariableIndices, boundSigns, boundValues,
                                           boundMode,
                                           maxNLIterations, tolr, tola,
                                           alpha, maxLSIterations, resetOuts );
    
                }
    
                if ( fatalErrorFlag ){
                    errorOut result = new errorNode( "homotopySolver", "Fatal error in Newton Raphson solution" );
                    result->addNext( error );
                    return result;
                }
                else if ( ( !convergeFlag ) && ( ds / 2 >= dsMin ) ){
                    s -= ds;
                    ds = std::max( ds / 2, dsMin );
                    s += ds;
    
                    floatOuts = oldFloatOuts;
                    intOuts   = oldIntOuts;
                }
                else if ( ( !convergeFlag ) && ( ds / 2 < dsMin ) ){
                    errorOut result = new errorNode( "homotopySolver", "Homotopy solver did not converge" );
                    result->addNext( error );
                    return result;
                }
            }

            xh = x;
        }

        //Solver completed successfully
        return NULL;
    }

    errorOut checkTolerance( const floatVector &R, const floatVector &tol, bool &result){
        /*!
         * Check whether the residual vector meets the tolerance returning a boolean.
         *
         * \param &R: The residual vector.
         * \param &tol: The tolerance.
         * \param result: The result
         */

        if (R.size() != tol.size()){
            return new errorNode("checkTolerance", "The residual and tolerance vectors don't have the same size");
        }

        for (unsigned int i=0; i<R.size(); i++){
            if (fabs(R[i])>tol[i]){
                result = false;
                return NULL;
            }
        }
        result = true;
        return NULL;
    }

    errorOut checkLSCriteria( const floatVector &R, const floatVector &Rp, bool &result, const floatType alpha){
        /*!
         * Perform the check on the line-search criteria setting result to false if the new residual does not meet it.
         * \f$l2norm(R) < (1 - alpha)* l2norm(Rp)\f$
         *
         * \param &R: The trial residual.
         * \param &Rp: the previous acceptable residual
         * \param &result: The output value.
         * \param &alpha: The scaling factor on Rp
         */

        if (R.size() != Rp.size()){
            return new errorNode("errorOut", "R and Rp have different sizes");
        }

        if (R.size() == 0){
            return new errorNode("errorOut", "R has a size of zero");
        }

        result = vectorTools::dot(R, R) < (1 - alpha)*vectorTools::dot(Rp, Rp);

        return NULL;
    }

    errorOut finiteDifference( stdFncNLF fxn, const floatVector &x0,
                            floatMatrix &grad, const floatMatrix &floatArgs, const intMatrix &intArgs, const floatType eps){
        /*!
         * Perform a forward finite difference gradient solve of the provided function.
         * Note that for functions that are more (or less) complex than this you may need to
         * wrap the function.
         *
         * The function function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to function
         * \param &intArgs: Additional integer arguments to the function
         * \param &value: The output value vector
         *
         * The main routine accepts the following parameters:
         * \param &x0: The initial iterate of x.
         * \param &grad: The finite difference gradient.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param eps: The perturbation. delta[i] = eps*(x0[i]) + eps
         */

        //Initialize the first value and the gradient
        floatVector y0, yi;

        //Compute the first value
        errorOut error = fxn(x0, floatArgs, intArgs, y0);
        if (error){
            errorOut result = new errorNode("finiteDifference", "Error in initial function calculation");
            result->addNext(error);
            return result;
        }

        grad = floatMatrix(y0.size(), floatVector(x0.size(), 0));

        for (unsigned int i=0; i<x0.size(); i++){
            //Set the step size
            floatVector delta = floatVector(x0.size(), 0);
            delta[i] = eps*fabs(x0[i]) + eps;

            //Compute the function after perturbation
            error = fxn(x0 + delta, floatArgs, intArgs, yi);
            if (error){
                errorOut result = new errorNode("finiteDifference", "Error in function calculation");
                result->addNext(error);
                return result;
            }

            //Set the terms of the gradient
            for (unsigned int j=0; j<yi.size(); j++){
                grad[j][i] = (yi[j] - y0[j])/delta[i];
            }
        }
        return NULL;
    }

    errorOut checkJacobian( stdFncNLFJ residual,
                            const floatVector &x0,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, bool &isGood, const floatType eps,
                            const floatType tolr, const floatType tola, const bool suppressOutput){
        /*!
         * Check if the jacobian is correct. Used as a debugging tool.
         *
         * The residual function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &residual: The residual vector
         * \param &jacobian: The jacobian matrix
         * \param &floatOuts: Additional returning floating point values.
         * \param &intOuts: Additional return integer values.
         *
         * The main routine accepts the following parameters:
         * \param &x0: The initial iterate of x.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param isGood: Whether the error in the jacobian is within tolerance.
         * \param eps: The perturbation. \f$delta[i] = eps*(x0[i]) + eps\f$
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param suppressOutput: Suppress the output to the terminal
         */

        //Wrap the residual function to hide the jacobian
        stdFncNLF residual_;
        residual_ = [&](const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                            floatVector &r){
            floatMatrix Jtmp;
            floatMatrix fO;
            intMatrix iO;
            return residual(x_, floatArgs_, intArgs_, r, Jtmp, fO, iO);
        };

        //Compute the finite difference jacobian
        floatMatrix finiteDifferenceJ;
        errorOut error = finiteDifference( residual_, x0, finiteDifferenceJ, floatArgs, intArgs);

        if (error){
            errorOut result = new errorNode("checkJacobian", "Error in finite difference");
            result->addNext(error);
            return error;
        }

        //Compute the analytic jacobian
        floatVector rtmp;
        floatMatrix analyticJ;
        floatMatrix floatOuts;
        intMatrix intOuts;
        residual(x0, floatArgs, intArgs, rtmp, analyticJ, floatOuts, intOuts);

        isGood = vectorTools::fuzzyEquals(finiteDifferenceJ, analyticJ, tolr, tola);

        if ((!isGood) && (!suppressOutput)){
            std::cout << "Jacobian is not within tolerance.\nError:\n";
            vectorTools::print(analyticJ - finiteDifferenceJ);
        }

        return NULL;
    }

    errorOut aFxn( const floatType &pseudoT, const floatType logAMax, floatType &a ){
        /*!
         * Compute the a parameter for the Barrier Function
         *
         * \param &pseudoT: The pseudo time ( 0 - 1 )
         * \param logAMax: The logarithm of the maximum a parameter value.
         * \param &a: The current value of a
         */

        a = std::exp( pseudoT * logAMax );

        return NULL;
    }

    errorOut aFxn( const floatType &pseudoT, const floatType logAMax, floatType &a, floatType &dadt ){
        /*!
         * Compute the a parameter for the Barrier Function along with the derivative w.r.t.
         * the pseudo time ( \f$t^s\f$ ).
         * 
         * \f$a = exp( log(A^{max}) t^s )
         *
         * \param &pseudoT: The pseudo time ( 0 - 1 )
         * \param logAMax: The logarithm of the maximum 'a' parameter value.
         * \param &a: The current value of 'a'
         * \param &dadt: The Jacobian of a w.r.t. pseudoT.
         */

        errorOut error = aFxn( pseudoT, logAMax, a );

        if ( error ){
            errorOut result = new errorNode( "aFxn (jacobian)", "Error in computation of a" );
            result->addNext( error );
            return result;
        }

        dadt = logAMax * a;

        return NULL;
    }

    errorOut computeBarrierFunction( const floatType &x, const floatType &pseudoT, const floatType &logAmax,
                                     const floatType &b, const bool &sign, floatType &barrierFunction ){
        /*!
         * Compute the barrier function for the constraint.
         *
         * \f$b = exp( s a( t^s ) ( b - x ) ) - 1\f$
         *
         * where
         * 
         * - \f$s =  1\f$ is a negative boundary
         * - \f$s = -1\f$ is a positive boundary
         *
         * \param &x: The constrained variable value ( \f$x\f$ ).
         * \param &pseudoT: The value of the pseudo time ( \f$t^s\f$ ).
         * \param &logAmax: The log of the maximum value of the \f$a\f$ parameter.
         * \param &b: The offset variable ( i.e. the location of the barrier ) ( \f$b\f$ )
         * \param &sign: A boolean which indicates if this is a positive ( 1 ) or
         *     negative ( 0 ) boundary ( \f$s\f$ ).
         * \param &barrierFunction: The value of the barrier function ( \f$b\f$ ).
         */

        floatType a;
        errorOut error = aFxn( pseudoT, logAmax, a );

        if ( error ){
            errorOut result = new errorNode( "computeBarrierFunction", "Error in the computation of the 'a' value" );
            result->addNext( error );
            return result;
        }

        floatType s = 1;
        if ( sign ){
            s = -1;
        }

        barrierFunction = std::exp( s * a * ( b - x ) ) - 1;

        return NULL;
    }

    errorOut computeBarrierFunction( const floatType &x, const floatType &pseudoT, const floatType &logAmax,
                                      const floatType &b, const bool &sign, floatType &barrierFunction,
                                      floatType &dbdx, floatType &dbdt ){
        /*!
         * Compute the barrier function for a positivity constraint where the barrier is defined as
         *
         * \f$b = exp( s * a * ( b - x ) ) - 1\f$
         *
         * where
         * 
         * - \f$s =  1\f$ is a negative boundary
         * - \f$s = -1\f$ is a positive boundary
         *
         * \param &x: The constrained variable value ( \f$x\f$ ).
         * \param &pseudoT: The value of the pseudo time ( \f$t^s\f$ ).
         * \param &logAmax: The log of the maximum value of the \f$a\f$ parameter.
         * \param &b: The offset variable ( i.e. the location of the barrier ) ( \f$b\f$ )
         * \param &sign: A boolean which indicates if this is a positive ( 1 ) or
         *     negative ( 0 ) boundary ( \f$s\f$ ).
         * \param &barrierFunction: The value of the barrier function ( \f$b\f$ ).
         * \param &dbdx: The Jacobian of the barrier function w.r.t. the variable value.
         * \param &dbdt: the Jacobian of the barrier function w.r.t. the pseudo time.
         */

        floatType a, dadt;
        errorOut error = aFxn( pseudoT, logAmax, a, dadt );

        if ( error ){
            errorOut result = new errorNode( "computeBarrierFunction (jacobian)", "Error in the computation of the 'a' value" );
            result->addNext( error );
            return result;
        }

        floatType s = 1;
        if ( sign ){
            s = -1;
        }

        barrierFunction = std::exp( s * a * ( b - x ) ) - 1;

        dbdx = -s * a * std::exp( s * a * ( b - x ) );
        dbdt = s * ( b - x ) * std::exp( s * a * ( b - x ) ) * dadt;

        return NULL;
    }

    errorOut computeBarrierHomotopyResidual( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                                                      floatVector &, floatMatrix &, floatMatrix &, intMatrix &
                                                                    ) > computeOriginalResidual,
                                             const solverTools::floatVector &x,
                                             const solverTools::floatMatrix &floatArgs, const solverTools::intMatrix &intArgs,
                                             solverTools::floatVector &residual, solverTools::floatMatrix &jacobian,
                                             solverTools::floatMatrix &floatOuts, solverTools::intMatrix &intOuts
                                           ){
        /*!
         * Compute the residual function for the barrier homotopy approach. This approach allows the user
         * to define barrier functions which enable a bounded root finding approach which can be very useful
         * when roots outside of the desired solution space have stronger basins of attraction than the
         * desired roots.
         * 
         * The updated residual is defined as
         * 
         * \f$R^b = \sum_{i=1}^{N^{barriers}} \left[\left( 1 - \frac{1}{a(t^s)} \right) R + \frac{1}{a\left(t^s\right)} b \right]\f$
         * 
         * where \f$R^b\f$ is the barrier residual, \f$N^{barriers}\f$ are the number of barrier functions to add,
         * \f$a\f$ is a parameter that is a function of pseudo-time \f$t^s\f$, and \f$b\f$ is the barrier function.
         *
         * \warning \b \emoji :warning: \emoji :warning: \emoji :warning: WARNING \emoji :warning: \emoji :warning: \emoji :warning:
         *     WARNING: If two `residualIndices` are identical, then only the second one will be used and the  first equation will not be observed.
         *
         * \todo{Add two-sided boundaries.}
         *
         * \param &x: The solution vector.
         * \param &floatArgs: The floating point arguments.
         * \param &intArgs: The integer arguments.
         * \param &residual: The residual vector.
         * \param &jacobian: The Jacobian matrix.
         * \param &floatOuts: The additional floating-point outputs.
         * \param &intOuts: The additional integer outputs.
         *
         * Note that floatArgs is modified such that
         * `floatArgs[ 0 ][ 0 ]`   = `pseudoTime` ( the homotopy pseudo-time )
         * `floatArgs[ 1 ]`        = `barrierValues` ( the values at which the barrier function activates )
         * `floatArgs[ 2 ]`        = `logAMaxValues` ( the maximum values of the \f$a\f$ parameter for the barrier function )
         * `floatArgs[ 1 -> end ]` = originalResidual `floatArgs` ( the `floatArgs` of the original residual function )
         *
         * Note that intArgs is modified such that
         * `intArgs[ 0 ]` = `variableIndices` ( the indices of the variables which have the barrier functions applied )
         * `intArgs[ 1 ]` = `residualIndices` ( the indices of the residual vector at which the barrier functions are applied )
         * `intArgs[ 2 ]` = `barrierSigns` ( the signs of the barrier functions 0 for negative barrier
         *     ( lower boundary ) and 1 for a positive barrier ( upper boundary )
         * `intArgs[ 3 -> end ]` = originalResidual `intArgs` ( the `intArgs` of the original residual function )
         *
         * The pseudo-time allows us to change the influence of the barrier function on the output.
         */

        if ( floatArgs.size() < 3 ){
            return new errorNode( "computeBarrierHomotopyResidual", "floatArgs must have at least a size of 3" );
        }

        if ( intArgs.size() < 3 ){
            return new errorNode( "computeBarrierHomotopyResidual", "intArgs must have at least a size of 3" );
        }

        //Extract the floatArgs values
        floatType pseudoTime      = floatArgs[ 0 ][ 0 ];
        floatVector barrierValues = floatArgs[ 1 ];
        floatVector logAMaxValues = floatArgs[ 2 ];

        floatMatrix floatArgsOriginalResidual( floatArgs.begin() + 3, floatArgs.begin() + floatArgs.size() );

        //Extract the intArgs values
        intVector variableIndices = intArgs[ 0 ];
        intVector residualIndices = intArgs[ 1 ];
        std::vector< bool > barrierSigns( intArgs[ 2 ].size() );
        for ( unsigned int i = 0; i < barrierSigns.size(); i++ ){
            barrierSigns[ i ] = ( bool )intArgs[ 2 ][ i ];
        }

        intMatrix intArgsOriginalResidual( intArgs.begin() + 3, intArgs.begin() + intArgs.size() );

        unsigned int nBarriers = variableIndices.size();

        if ( ( residualIndices.size() != nBarriers ) || ( barrierValues.size() != nBarriers ) || ( logAMaxValues.size() != nBarriers ) || ( barrierSigns.size() != nBarriers ) ){
            std::string output_message = "The sizes of variableIndices, residualIndices, barrierValues, and logAMaxValues are not the same\n";
            output_message            += "    variableIndices: " + std::to_string( variableIndices.size() ) + "\n";
            output_message            += "    residualIndices: " + std::to_string( residualIndices.size() ) + "\n";
            output_message            += "    barrierValues:   " + std::to_string( barrierValues.size() ) + "\n";
            output_message            += "    logAMaxValues:   " + std::to_string( logAMaxValues.size() ) + "\n";
            output_message            += "    barrierSigns:    " + std::to_string( barrierSigns.size() ) + "\n";
            return new errorNode( "computeBarrierHomotopyResidual", output_message.c_str() );
        }

        //Evaluate the original residual
        floatVector originalResidual;
        floatMatrix originalJacobian;
        floatMatrix originalFloatOuts = floatOuts;
        errorOut error = computeOriginalResidual( x,
                                                  floatArgsOriginalResidual, intArgsOriginalResidual,
                                                  originalResidual, originalJacobian,
                                                  originalFloatOuts, intOuts
                                                );
        residual = originalResidual;
        jacobian = originalJacobian;

        if ( error ){
            errorOut result = new errorNode( "computeBarrierHomotopyResidual", "Error in the computation of the plastic deformation residual" );
            result->addNext( error );
            return result;
        }

        //Compute the values of the barrier functions and the weighting functions
        floatType barrierFunction = 0;
        floatType dbdx = 0;
        floatType dbdt = 0;
        floatType a    = 1;
        floatType dadt = 0;

        //Update the gradients of the residuals and Jacobian
        floatVector dresidualdt( residual.size(), 0 );

        for ( unsigned int i = 0; i < nBarriers; i++ ){
            error = computeBarrierFunction( x[ variableIndices[ i ] ], pseudoTime, logAMaxValues[ i ],
                                            barrierValues[ i ], barrierSigns[ i ],
                                            barrierFunction, dbdx, dbdt );

            if ( error ){
                std::string output_message = "Error in the computation of barrier function " + std::to_string( i );
                errorOut result = new errorNode( "computeBarrierHomotopyResidual", output_message.c_str() );
                result->addNext( error );
                return result;
            }

            //Compute the weighting values
            error = aFxn( pseudoTime, logAMaxValues[ i ], a, dadt );

            if ( error ){
                std::string output_message = "Error in the computation of the 'a' parameter of barrier equation " + std::to_string( i );
                errorOut result = new errorNode( "computeBarrierHomotopyResidual", output_message.c_str() );
                result->addNext( error );
                return result;
            }

            //Assemble the homotopy residual
            residual[ residualIndices[ i ] ] = ( 1. - 1. / a ) * originalResidual[ residualIndices[ i ] ] + ( 1. / a ) * barrierFunction;

            //Assemble the derivative of the homotopy residual w.r.t. the pseudo time
            dresidualdt[ residualIndices[ i ] ] = 1 / ( a * a ) * dadt * originalResidual[ residualIndices[ i ] ]
                                                - 1 / ( a * a ) * dadt * barrierFunction
                                                + ( 1 / a ) * dbdt;

            //Add the terms to the jacobian ( dresidualdx )
            for ( unsigned int j = 0; j < jacobian[ i ].size(); j++ ){
                jacobian[ residualIndices[ i ] ][ j ] = ( 1. - 1. / a ) * originalJacobian[ residualIndices[ i ] ][ j ];
            }

            jacobian[ residualIndices[ i ] ][ variableIndices[ i ] ] += ( 1. / a ) * dbdx;
        }

        //Save the jacobian of the residual w.r.t. the pseudo time. This is inserted at the beginning of the floatOuts.
        floatOuts = floatMatrix( originalFloatOuts.size() + 1 );
        floatOuts[ 0 ] = dresidualdt;
        for ( unsigned int i = 0; i < originalFloatOuts.size(); i++ ){
            floatOuts[ i + 1 ] = originalFloatOuts[ i ];
        }

        return NULL;
    }

    errorOut applyBoundaryLimitation( const floatVector &x0, const intVector &variableIndices, const intVector &barrierSigns,
                                      const floatVector &barrierValues, floatVector &dx, const floatType tolr,
                                      const floatType tola, const bool mode ){
        /*!
         * Apply the boundary limitation to the update step.
         * dx will either be scaled so that all of the variables respect the boundaries or each variable which
         * violates the constraint will be set to the boundary value. This is should be viewed as
         * a last resort if even homotopy has failed as it is increasingly likely that it will
         * not result in convergence.
         *
         * \param &x0: The base vector from which dx extends.
         * \param &variableIndices: The indices of the vector that have constraints applied to them.
         * \param &barrierSigns: The sign of the barrier 0 is a negative barrier ( i.e. lower bound )
         *     1 is a positive barrier ( i.e. upper bound )
         * \param &barrierValues: The locations of the barriers.
         * \param &dx: The change in x vector.
         * \param &tolr: The relative tolerance. Defaults to 1e-9
         * \param &tola: The absolute tolerence. Defautls to 1e-9
         * \param mode: The mode of the boundary limitation. If false, all of dx is scaled, if true
         *     only the value that violates the constraint is set to the constraint.
         */

        unsigned int nBounds = variableIndices.size();

        if ( ( barrierSigns.size() != nBounds ) || ( barrierValues.size() != nBounds ) ){

            std::string output_message = "The defined barrier are not consistent in size\n";
            output_message            += "    variableIndices: " + std::to_string( variableIndices.size() ) + "\n";
            output_message            += "    barrierSigns:    " + std::to_string( barrierSigns.size() ) + "\n";                                                                                                            output_message            += "    barrierValues:   " + std::to_string( barrierValues.size() ) + "\n";
            return new errorNode( "applyBoundaryLimitations", output_message.c_str() );
        }

        //Initialize variables
        int index;
        floatType tol, d;
        floatType scaleFactor = 1.0;

        //Do the barrier search
        for ( unsigned int i = 0; i < nBounds; i++ ){

            //Extract the index of the variable
            index = variableIndices[ i ];
 
            //Calculate the tolerance
            //
            //This relative tolerance is not the only tolerance that could be used.
            //It does allow for the constraints to be violated but this should be very small
            //compared to the magnitude of the previously converged value.                              
            tol = tolr * fabs( x0[ index ] ) + tola;

            //Determine the amount of constraint violation
            d = ( x0[ index ] + dx[ index ] ) - barrierValues[ i ];
            if ( barrierSigns[ i ] == 0 ){
                d *= -1;
            }
            else if ( barrierSigns[ i ] != 1 ){
                return new errorNode( "applyBoundaryLimitation", "The barrier sign must be zero or 1" );
            }

            d = 0.5 * ( d + fabs( d ) );

            //Determine the required scale factor
            if ( d > tol ){

                if ( mode ){

                    dx[ index ] = barrierValues[ i ] - x0[ index ];

                }
                else{

                    scaleFactor = std::min( scaleFactor, ( barrierValues[ i ] - x0[ index ] ) / dx[ index ] );
                }
            }
        }

        //Scale the dx vector if appropriate
        if ( !mode ){
            if ( !vectorTools::fuzzyEquals( scaleFactor, 1.0 ) ){
                dx *= scaleFactor;
            }
        }

        return NULL;
    }

    errorOut barrierHomotopySolver( stdFncNLFJ residual, const floatType &dt, const floatVector &x0,
                                    const intVector &variableIndices, const intVector &residualIndices,
                                    const intVector &barrierSigns, const floatVector &barrierValues,
                                    const floatVector &logAMaxValues,
                                    const floatMatrix &floatArgs, const intMatrix &intArgs,
                                    const bool &implicitRefine,
                                    floatVector &x, bool &convergeFlag, bool &fatalErrorFlag,
                                    floatMatrix &floatOuts, intMatrix &intOuts,
                                    const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                                    const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * Perform a non-linear solve using a homotopy method with barriers.
         *
         * Note that if barriers are added to a variable multiple times or to a single component of the residual
         * equation several times unexpected responses can result.
         *
         * \param residual: The residual function to initialize
         * \param &dt: The pseudo-time step size.
         * \param &x0: The initial iterate of the solution vector. It is strongly suggested to
         *     use a vector which
         * \param &variableIndices: The indices of the variables which are to have barrier conditions
         *     applied.
         * \param &residualIndices: The indices of the residual vector which should have the barrier
         *     equations applied.
         * \param &barrierSigns: The signs of the barriers.
         * \param &barrierValues: The locations of the barriers.
         * \param &logAMaxValues: The log of the maximum values of the a parameter.
         *     these should be as large as possible without causing numeric issues. Values of 10 seem to work well
         *     for the exponential barrier.
         * \param &floatArgs: The floating point arguments for the residual equation.
         * \param &intArgs: The integar arguments for the residual equation.
         * \param &implicitRefine: Boolean which indicates if an implicit refining of the
         *     explicit pseudo-time step should occur. This can add significant computational expense
         *     and should only be used if necessary.
         * \param &x0: The initial iterate of the solution vector.
         * \param &x: The solution vector.
         * \param &convergeFlag: The flag which indicates convergence.
         * \param &fatalErrorFlag: The flag which indicates the presence of fatal errors.
         * \param &floatOuts: The floating point outputs for the residual equation.
         * \param &intOuts: The integer outputs for the residual equation.
         * \param maxNLIterations: The maximum number of non-linear iterations for the Newton-Raphson
         *     solve.
         * \param tolr: The relative tolerance for the Newton-Raphson solve.
         * \param tola: The absolute tolerance for the Newton-Raphson solve.
         * \param alpha: The alpha parameter for the line search.
         * \param maxLSIterations: The maximum number of line search iterations.
         * \param resetOuts: The flag for whether the output matrices should be reset
         *     prior to each iteration.
         */

        solverType linearSolver;
        floatMatrix jacobian;
        return barrierHomotopySolver( residual, dt, x0, variableIndices, residualIndices, barrierSigns, barrierValues,
                                      logAMaxValues, floatArgs, intArgs, implicitRefine, x, convergeFlag, fatalErrorFlag,
                                      floatOuts, intOuts, linearSolver, jacobian,
                                      maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts );
    }

    errorOut barrierHomotopySolver( stdFncNLFJ residual, const floatType &dt, const floatVector &x0,
                                    const intVector &variableIndices, const intVector &residualIndices,
                                    const intVector &barrierSigns, const floatVector &barrierValues,
                                    const floatVector &logAMaxValues,
                                    const floatMatrix &floatArgs, const intMatrix &intArgs,
                                    const bool &implicitRefine,
                                    floatVector &x, bool &convergeFlag, bool &fatalErrorFlag,
                                    floatMatrix &floatOuts, intMatrix &intOuts, solverType &linearSolver, floatMatrix &jacobian,
                                    const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                                    const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * Perform a non-linear solve using a homotopy method with barriers.
         *
         * Note that if barriers are added to a variable multiple times or to a single component of the residual
         * equation several times unexpected responses can result.
         *
         * \param residual: The residual function to initialize
         * \param &dt: The pseudo-time step size.
         * \param &x0: The initial iterate of the solution vector. It is strongly suggested to
         *     use a vector which
         * \param &variableIndices: The indices of the variables which are to have barrier conditions
         *     applied.
         * \param &residualIndices: The indices of the residual vector which should have the barrier
         *     equations applied.
         * \param &barrierSigns: The signs of the barriers.
         * \param &barrierValues: The locations of the barriers.
         * \param &logAMaxValues: The log of the maximum values of the a parameter.
         *     these should be as large as possible without causing numeric issues. Values of 10 seem to work well
         *     for the exponential barrier.
         * \param &floatArgs: The floating point arguments for the residual equation.
         * \param &intArgs: The integar arguments for the residual equation.
         * \param &implicitRefine: Boolean which indicates if an implicit refining of the
         *     explicit pseudo-time step should occur. This can add significant computational expense
         *     and should only be used if necessary.
         * \param &x0: The initial iterate of the solution vector.
         * \param &x: The solution vector.
         * \param &convergeFlag: The flag which indicates convergence.
         * \param &fatalErrorFlag: The flag which indicates the presence of fatal errors.
         * \param &floatOuts: The floating point outputs for the residual equation.
         * \param &intOuts: The integer outputs for the residual equation.
         * \param &linearSolver: The linear solver object.
         * \param &J: The Jacobian matrix.
         * \param maxNLIterations: The maximum number of non-linear iterations for the Newton-Raphson
         *     solve.
         * \param tolr: The relative tolerance for the Newton-Raphson solve.
         * \param tola: The absolute tolerance for the Newton-Raphson solve.
         * \param alpha: The alpha parameter for the line search.
         * \param maxLSIterations: The maximum number of line search iterations.
         * \param resetOuts: The flag for whether the output matrices should be reset
         *     prior to each iteration.
         */

        //Initialize the pseudo time
        floatType pseudoTime = 0.;

        //Initialize the output error
        errorOut error;

        //Form argument matrices
        floatMatrix homotopyFloatArgs( floatArgs.size() + 3 );
        intMatrix   homotopyIntArgs( intArgs.size() + 3 );

        homotopyFloatArgs[ 0 ] = { pseudoTime };
        homotopyFloatArgs[ 1 ] = barrierValues;
        homotopyFloatArgs[ 2 ] = logAMaxValues;

        for ( unsigned int i = 0; i < floatArgs.size(); i++ ){
            homotopyFloatArgs[ 3 + i ] = floatArgs[ i ];
        }

        homotopyIntArgs[ 0 ] = variableIndices;
        homotopyIntArgs[ 1 ] = residualIndices;
        homotopyIntArgs[ 2 ] = barrierSigns;

        for ( unsigned int i = 0; i < intArgs.size(); i++ ){
            homotopyIntArgs[ 3 + i ] = intArgs[ i ];
        }

        //Form output matrices
        floatMatrix homotopyFloatOuts = floatOuts;
        intMatrix homotopyIntOuts = intOuts;

        //Wrap the barrier homotopy function
        stdFncNLFJ homotopyResidual;
        homotopyResidual = [&](const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                               floatVector &r, floatMatrix &J, floatMatrix &fO, intMatrix &iO ){

            error = computeBarrierHomotopyResidual( residual, x_, floatArgs_, intArgs_, r, J, fO, iO );

            if ( error ){
                errorOut result = new errorNode( "barrierHomotopySolver::homotopyResidual", "Error in wrapped barrier homotopy residual function" );
                result->addNext( error );
                return result;
            }

            return static_cast<errorOut>(NULL);
        };

        //Perform solve to initialize x0

        error = newtonRaphson( homotopyResidual, x0, x, convergeFlag, fatalErrorFlag, homotopyFloatOuts, homotopyIntOuts,
                               homotopyFloatArgs, homotopyIntArgs,
                               maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts
                             );

        if ( error ){
            errorOut result = new errorNode( "barrierHomotopySolver", "Error in initial Newton-Raphson solve" );
            result->addNext( error );
            return result;
        }

        //Initialize variables
        floatVector R, x0_update;
        floatMatrix J;
        unsigned int rank;

        while ( pseudoTime < 1.0 ){

            //Reset the additional outputs
            homotopyFloatOuts = floatOuts;
            homotopyIntOuts = intOuts;

            //Update the pseudo time
            homotopyFloatArgs[ 0 ][ 0 ] = pseudoTime;

            //Evaluate the homotopy residual
            error = computeBarrierHomotopyResidual( residual, x, homotopyFloatArgs, homotopyIntArgs,
                                                    R, J, homotopyFloatOuts, homotopyIntOuts );

            if ( error ){
                errorOut result = new errorNode( "barrierHomotopySolver", "Error in computation of homotopy residual" );
                result->addNext( error );
                return result;
            }

            //Solve for the new dt
            floatType dPT = dt;
            if ( pseudoTime + dt > 1. ){
                dPT = ( 1 - pseudoTime );
            }

            //Solve for the new dx
            x -= dPT * vectorTools::solveLinearSystem( J, homotopyFloatOuts[ 0 ], rank );

            //Update the pseudo-time
            pseudoTime += dt;
            if ( pseudoTime > 1.0 ){
                pseudoTime = 1.0;
            }

            //Refine answer if required. This runs an implicit solver to solve the
            //homotopy function which ensures that the derivatives used for the
            //next pseudo-time step are valid
            if ( implicitRefine ){
                homotopyFloatOuts = floatOuts;
                homotopyIntOuts = intOuts;

                homotopyFloatArgs[ 0 ][ 0 ] = pseudoTime;
                x0_update = x;

                error = newtonRaphson( homotopyResidual, x0_update, x, convergeFlag, fatalErrorFlag,
                                       homotopyFloatOuts, homotopyIntOuts,
                                       homotopyFloatArgs, homotopyIntArgs,
                                       maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts
                                     );

                if ( error ){
                    errorOut result = new errorNode( "barrierHomotopySolver", "Error in Newton-Raphson solve during implicit update" );
                    result->addNext( error );
                    return result;
                }
            }
        }

        //Using the initialized value of x, run a Newton-Raphson solver
        x0_update = x;
        error = newtonRaphson( residual, x0_update, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts,
                               floatArgs, intArgs, linearSolver, jacobian,
                               maxNLIterations, tolr, tola, alpha, maxLSIterations, resetOuts );

        if ( error ){
            errorOut result = new errorNode( "barrierHomotopySolver",
                                             "Error in the final Newton-Raphson solution of the barrier homotopy solver" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut BFGS( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                            floatType &, floatVector &, floatMatrix &, intMatrix &
                                          ) > lagrangianGradientFunction,
                   const floatVector &x0,
                   floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                   const floatMatrix &floatArgs, const intMatrix &intArgs,
                   const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                   const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts,
                   const floatType stepSize, const floatType maxdx ){
        /*!
         * An implementation of the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm for solving optimization problems.
         *
         * \param lagrangianFunction: The Lagrangian function
         * \param lagrangianGradientFunction: The gradient of the Lagrangian function.
         * \param &x0: The intial iterate of x
         * \param &x: The converged value of x
         * \param &convergeFlag: A flag indicating if convergence has been achieved.
         * \param &fatalErrorFlag: A flag indicating if a fatal error was encountered.
         * \param &floatOuts: A matrix of floating-point outputs.
         * \param &intOuts: A matrix of integer outputs.
         * \param &floatArgs: A matrix of floating-point arguments.
         * \param &intArgs: A matrix of integer arguments.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance.
         * \param tola: The absolute tolerance.
         * \param alpha: The line-search parameter.
         * \param maxLSIterations: The maximum number of line-search reductions.
         * \param resetOuts: The flag to indicate if the floating-point outputs should be reset at each
         *     iteration.
         * \param stepSize: The estimated step size. Defaults to 1.
         * \param maxdx: The maximum allowable step size in terms of the norm of the solution
         *     vector. If negative this is ignored.
         *
         * The lagrangianGradient function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param intMatrix &intArgs: Additional integer arguments to the residual
         * \param &value: The value of the Lagrangian.
         * \param &gradient: The gradient of the Lagrangian function.
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         */

        //Solve for the initial gradient of the Lagrangian
        floatType lagrangian_k, lagrangian_kp1;
        floatVector lagrangianGradient_k, lagrangianGradient_kp1;

        floatMatrix floatOuts0;
        intMatrix intOuts0;

        if ( resetOuts ){
            floatOuts0 = floatOuts;
            intOuts0 = intOuts;
        }

        errorOut error = lagrangianGradientFunction( x0, floatArgs, intArgs, lagrangian_kp1, lagrangianGradient_kp1,
                                                     floatOuts, intOuts );

        if ( error ){
            errorOut result = new errorNode( "BFGS", "Error in computation of the Lagrangian Gradient" );
            result->addNext( error );
            convergeFlag = false;
            fatalErrorFlag = true;
            return result;
        }

        //Set the initial iterate of the hessian
        floatMatrix B = ( vectorTools::l2norm( lagrangianGradient_kp1 ) / stepSize ) * vectorTools::eye< floatType >( lagrangianGradient_kp1.size() );

        //Set the tolerance for each value individually
        floatVector tol = floatVector( lagrangianGradient_kp1.size(), 0 );
        for ( unsigned int i = 0; i < lagrangianGradient_kp1.size(); i++ ){
            tol[ i ] = tolr * fabs( lagrangianGradient_kp1[ i ] ) + tola;
        }

        //Check if convergence has been achieved
        bool converged;
        error = checkTolerance( lagrangianGradient_kp1, tol, converged );

        if ( error ){
            errorOut result = new errorNode( "BFGS", "Error in tolerence check" );
            result->addNext( error );
            fatalErrorFlag = true;
            return result;
        }

        //Copy lagrangianGradient_kp1 to lagrangianGradient_p
        lagrangian_k = lagrangian_kp1;
        lagrangianGradient_k = lagrangianGradient_kp1;

        //Initialize the iteration
        floatType theta;
        floatType pNorm;
        floatVector r, p, s, y, Bs;
        floatType ys;
        floatType sBs;

        floatMatrix floatOuts_k;
        intMatrix intOuts_k;

        x = x0;

        unsigned int niter = 0;
        unsigned int nLSIterations;
        //Begin the iteration loop
        while ( ( !converged ) && ( niter < maxNLIterations ) ){

            //Compute the direction
            unsigned int rank;
            p = -vectorTools::solveLinearSystem( B, lagrangianGradient_k, rank );
            pNorm = vectorTools::l2norm( p );

            //Apply maximum dx step size limitation
            if ( ( maxdx > 0 ) && ( pNorm > maxdx ) ){
                p *= ( maxdx / pNorm );
            }

            if ( rank < lagrangianGradient_k.size() ){
                convergeFlag = false;
                fatalErrorFlag = false;
                return new errorNode( "BFGS", "The approximate Hessian is singular" );
            }

            if ( resetOuts ){
                floatOuts = floatOuts0;
                intOuts = intOuts0;
            }

            //Save the current value of the outputs
            if ( resetOuts ){
                floatOuts_k = floatOuts0;
                intOuts_k = intOuts0;
                floatOuts = floatOuts0;
                intOuts = intOuts0;

            }
            else{
                floatOuts_k = floatOuts;
                intOuts_k = intOuts;
            }

            error = lagrangianGradientFunction( x + p, floatArgs, intArgs, lagrangian_kp1, lagrangianGradient_kp1,
                                                floatOuts, intOuts );

            if ( error ){
                errorOut result = new errorNode( "BFGS", "Error in computation of the Lagrangian gradient function" );
                result->addNext( error );
                fatalErrorFlag = true;
                return result;
            }

            //Begin the line search
            floatType lambda = 1;
            nLSIterations = 0;

            while ( ( lagrangian_kp1 > ( 1 - alpha ) * lagrangian_k ) && ( nLSIterations < maxLSIterations ) ){

                //Reduce the LS step-size
                lambda *= 0.5;

                floatOuts = floatOuts_k;
                intOuts = intOuts_k;

                error = lagrangianGradientFunction( x + lambda * p, floatArgs, intArgs, lagrangian_kp1, lagrangianGradient_kp1,
                                                    floatOuts, intOuts );

                if ( error ){
                    errorOut result = new errorNode( "BFGS", "Error in the computation of the Lagrangian gradient function of the line search" );
                    result->addNext( error );
                    fatalErrorFlag = true;
                    return result;
                }

                nLSIterations++;
            }

            if ( lagrangian_kp1 > ( 1 - alpha ) * lagrangian_k ){
                convergeFlag = false;
                fatalErrorFlag = false;
                return new errorNode( "BFGS", "Line search did not converge" );
            }
            else{

                //Continue with BFGS update

                //Compute s
                s = lambda * p;

                //Update x
                x += s;

                //Compute y
                y = lagrangianGradient_kp1 - lagrangianGradient_k;

                //Compute the damping factor
                theta = 1.;

                //Update the approximation of the hessian
                ys = vectorTools::dot( y, s );
                Bs = vectorTools::dot( B, s );
                sBs = vectorTools::dot( s, Bs );

                if ( ys < 0.2 * sBs ){
                    theta = 0.8 * sBs / ( sBs - ys );
                }

                r = theta * y  + ( 1 - theta ) * Bs;

                B += vectorTools::dyadic( r, r ) / vectorTools::dot( s, r )
                   - vectorTools::dyadic( Bs, Bs ) / vectorTools::dot( s, Bs );

                //Update the previous values
                lagrangian_k = lagrangian_kp1;
                lagrangianGradient_k = lagrangianGradient_kp1;

                //Check the convergence tolerence
                error = checkTolerance( lagrangianGradient_kp1, tol, converged );

                if ( error ){
                    return new errorNode( "BFGS", "Error in the tolerence check at end of iteration" );
                }
            }
        }

        if ( converged ){
            convergeFlag = true;
            fatalErrorFlag = false;
            return NULL;
        }
        else{
            convergeFlag = false;
            fatalErrorFlag = false;

            return new errorNode( "BFGS", "BFGS did not converge" );
        }
    }

    errorOut homotopyBFGS( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatType &, floatVector &, floatMatrix &, intMatrix &
                                                   ) > lagrangianGradientFunction,
                           const floatVector &x0,
                           floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                           const floatMatrix &floatArgs, const intMatrix &intArgs,
                           const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                           const floatType alpha, const unsigned int maxLSIterations, const floatType ds0,
                           const floatType dsMin, const bool resetOuts, const floatType maxdx ){
        /*!
         * Optimize a non-linear equation using a homotopy BFGS method. This method
         * can be successful in solving very stiff equations which other techniques
         * struggle to capture. It effectively breaks the solve up into sub-steps
         * of easier to solve equations which will eventually converge to the
         * more difficult problem.
         *
         * \warning \b \emoji :warning: \emoji :warning: \emoji :warning: WARNING \emoji :warning: \emoji :warning: \emoji :warning:
         *     WARNING: This function is less tested than would be desired and should be used with caution.
         *
         * The lagrangian function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &lagrangian: The residual vector
         * \param &lagrangianGradient: The jacobian matrix of the residual w.r.t. x
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         *
         * The main routine accepts the following parameters:
         * \param &x0: The initial iterate of x.
         * \param &x: The converged value of the solver.
         * \param &convergeFlag: A flag which indicates whether the solver converged.
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         * \param &floatArgs: The additional floating-point arguments.
         * \param &intArgs: The additional integer arguments.
         * \param maxNLIterations: The maximum number of non-linear iterations.
         * \param tolr: The relative tolerance
         * \param tola: The absolute tolerance
         * \param alpha: The line search criteria.
         * \param maxLSIterations: The maximum number of line-search iterations.
         * \param ds0: The initial pseudo-time step. Defaults to 1.0
         * \param dsMin: The minimum pseudo-time step size. Defaults to 0.1
         * \param resetOuts: The flag for whether the output matrices should be reset
         *     prior to each iteration.
         * \param const floatType maxdx: The maximum allowable change in the solution vector. If negative this is ignored.
         */

        //Initialize the homotopy solver
        floatType ds = ds0;
        floatType s  = 0;
        floatVector xh = x0;

        //Save the floatOuts and intOuts
        floatMatrix oldFloatOuts = floatOuts;
        intMatrix   oldIntOuts   = intOuts;

        //Initialize the error output
        errorOut error;


        //Define the homotopy residual
        stdFncLagrangianG homotopyLagrangianGradient;
        homotopyLagrangianGradient = [&]( const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                                          floatType &l, floatVector &dldx, floatMatrix &fO, intMatrix &iO ){

            floatType L;
            floatVector dLdx;

            error = lagrangianGradientFunction( x_, floatArgs_, intArgs_, L, dLdx, fO, iO );

            if (error){
                errorOut result = new errorNode("homotopyBFGS::homotopyLagrangianGradient", "error in lagrangian gradient calculation");
                result->addNext(error);
                return result;
            }

            l = 0.5 * ( 1 - s ) * vectorTools::dot( x_ - x0, x_ - x0 ) + s * L;
            dldx = ( 1 - s ) * ( x_ - x0 ) + s * dLdx;

            return static_cast<errorOut>(NULL);
        };

        //Begin the homotopy loop
        while ( s < 1 ){
            //Update s
            s += ds;
            s = std::min( s, 1. );

            //Initialize the solver
            convergeFlag = false;

            if ( !resetOuts ){
                oldFloatOuts = floatOuts;
                oldIntOuts   = intOuts;
            }
            else{
                floatOuts = oldFloatOuts;
                intOuts = oldIntOuts;
            }

            //Begin the adaptive homotopy loop
            while ( !convergeFlag ){

                error = BFGS( homotopyLagrangianGradient, xh, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts,
                              floatArgs, intArgs,
                              maxNLIterations, tolr, tola,
                              alpha, maxLSIterations, resetOuts, maxdx );

                if ( fatalErrorFlag ){
                    errorOut result = new errorNode( "homotopyBFGS", "Fatal error in Newton Raphson solution" );
                    result->addNext( error );
                    return result;
                }

                else if ( ( !convergeFlag ) && ( ds / 2 > dsMin ) ){
                    s -= ds;
                    ds = std::max( ds / 2, dsMin );
                    s += ds;

                    floatOuts = oldFloatOuts;
                    intOuts   = oldIntOuts;
                }

                else if ( ( !convergeFlag ) && ( ds / 2 < dsMin ) ){
                    errorOut result = new errorNode( "homotopyBFGS", "Homotopy solver did not converge" );
                    result->addNext( error );
                    return result;
                }
            }

            xh = x;
        }

        //Solver completed successfully
        return NULL;
    }

}
