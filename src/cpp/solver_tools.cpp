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

    errorOut newtonRaphson( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                            const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation
         * of a typical Newton-Raphson solver which can take an arbitrary
         * residual function.
         *
         * The residual function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &residual: The residual vector
         * \param &jacobian: The jacobian matrix of the residual w.r.t. x
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         *
         * The main routine accepts the following parameters:
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

    errorOut newtonRaphson( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                            const floatVector &x0,
                            floatVector &x, bool &convergeFlag, bool &fatalErrorFlag, floatMatrix &floatOuts, intMatrix &intOuts,
                            const floatMatrix &floatArgs, const intMatrix &intArgs, solverType &linearSolver, floatMatrix &J,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const bool resetOuts ){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation
         * of a typical Newton-Raphson solver which can take an arbitrary
         * residual function.
         *
         * The residual function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &residual: The residual vector
         * \param &jacobian: The jacobian matrix of the residual w.r.t. x
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         *
         * The main routine accepts the following parameters:
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

    errorOut newtonRaphson( std::function< errorOut( const floatVector &, const floatMatrix &, const intMatrix &,
                                                     floatVector &, floatMatrix &, floatMatrix &, intMatrix & ) > residual,
                            const floatVector &x0,
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
         * The residual function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &residual: The residual vector
         * \param &jacobian: The jacobian matrix of the residual w.r.t. x
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         *
         * The main routine accepts the following parameters:
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
         *     more details..
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

    errorOut homotopySolver( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                             const floatVector &x0,
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
         * The residual function should have inputs of the form
         * \param &x: A vector of the variable to be solved.
         * \param &floatArgs: Additional floating point arguments to residual
         * \param &intArgs: Additional integer arguments to the residual
         * \param &residual: The residual vector
         * \param &jacobian: The jacobian matrix of the residual w.r.t. x
         * \param &floatOuts: Additional floating point values to return.
         * \param &intOuts: Additional integer values to return.
         *
         * The main routine accepts the following parameters:
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

        //Initialize the homotopy solver
        floatType ds = ds0;
        floatType s  = 0;
        floatVector xh = x0;
        floatVector xdot = x0;
        floatVector rh, dxh;
        solverType ls;
        unsigned int rank;

        //Save the flowOuts and intOuts
        floatMatrix oldFloatOuts = floatOuts;
        intMatrix   oldIntOuts   = intOuts;

        //Define the homotopy residual equation
        floatVector Rinit;
        floatMatrix J;

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
        homotopyResidual = [&]( const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                                floatVector &r, floatMatrix &_J, floatMatrix &fO, intMatrix &iO ){

            floatVector R;

            error = residual( x_, floatArgs_, intArgs_, R, _J, fO, iO );

            if ( error ){
                r = R;
                errorOut result = new errorNode( "homotopySolver::homotopyResidual", "error in residual calculation" );
                result->addNext( error );
                return result;
            }

            r = R - ( 1 - s ) * Rinit;

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

                //Compute the eplicit estimate of xh ( this is kind of what makes it a homotopy )
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
    
                    error = newtonRaphson( homotopyResidual, xh + dxh, x, convergeFlag, fatalErrorFlag, floatOuts, intOuts,
                                           floatArgs, intArgs, maxNLIterations, tolr, tola,
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
         * /param &x0: The base vector from which dx extends.
         * /param &variableIndices: The indices of the vector that have constraints applied to them.
         * /param &barrierSigns: The sign of the barrier 0 is a negative barrier ( i.e. lower bound )
         * /   1 is a positive barrier ( i.e. upper bound )
         * /param &barrierValues: The locations of the barriers.
         * /param &dx: The change in x vector.
         * /param &tolr: The relative tolerance. Defaults to 1e-9
         * /param &tola: The absolute tolerence. Defautls to 1e-9
         * /param mode: The mode of the boundary limitation. If false, all of dx is scaled, if true
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

}
