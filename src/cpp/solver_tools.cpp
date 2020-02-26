/*!
===============================================================================
|                             solver_tools.cpp                                |
===============================================================================
| The solver tools library. Incorporates a collection of non-linear solver    |
| tools built on top of Eigen. These tools are intended to be general enough  |
| to solve any small-ish nonlinear problem that can be solved using a Newton  |
| Raphson approach.                                                           |
===============================================================================
*/

#include<solver_tools.h>

namespace solverTools{

    errorOut newtonRaphson( std::function< errorOut(const floatVector &, const floatMatrix &, const intMatrix &,
                                                    floatVector &, floatMatrix &, floatMatrix &, intMatrix &) > residual,
                            const floatVector &x0,
                            floatVector &x, bool &convergeFlag, floatMatrix &floatOuts, intMatrix &intOuts, 
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations){
        /*!
         * The main Newton-Raphson non-linear solver routine. An implementation 
         * of a typical Newton-Raphson solver which can take an arbitrary 
         * residual function.
         * 
         * The residual function should have inputs of the form
         * :param const floatVector &x: A vector of the variable to be solved.
         * :param const floatMatrix &floatArgs: Additional floating point arguments to residual
         * :param const intMatrix &intArgs: Additional integer arguments to the residual
         * :param floatVector &residual: The residual vector
         * :param floatMatrix &jacobian: The jacobian matrix of the residual w.r.t. x
         * :param floatMatrix &floatOuts: Additional floating point values to return.
         * :param intMatrix &intOuts: Additional integer values to return.
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param floatVector &x: The converged value of the solver.
         * :param bool &convergeFlag: A flag which indicates whether the solver converged.
         * :param floatMatrix &floatOuts: Additional floating point values to return.
         * :param intMatrix &intOuts: Additional integer values to return.
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         * :param const unsigned int maxNLIterations: The maximum number of non-linear iterations.
         * :param const floatType tolr: The relative tolerance
         * :param const floatType tola: The absolute tolerance 
         * :param const floatType alpha: The line search criteria.
         * :param const unsigned int maxLSIterations: The maximum number of line-search iterations.
         */

        //Compute the initial residual and jacobian
        floatVector dx = floatVector(x0.size(), 0);
        floatVector ddx;
        floatVector R, Rp;
        floatMatrix J;

        errorOut error = residual(x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts);

        if (error){
            errorOut result = new errorNode("newtonRaphson", "Error in computation of initial residual");
            result->addNext(error);
            return result;
        }

        if (R.size() != x0.size()){
            return new errorNode("newtonRaphson", "The residual and x0 don't have the same lengths. The problem is ill-defined.");
        }

        //Set the tolerance for each value individually
        floatVector tol = floatVector(R.size(), 0);
        for (unsigned int i=0; i<R.size(); i++){tol[i] = tolr*fabs(R[i]) + tola;}

        //Copy R to Rp
        Rp = R;

        //Initialize variables required for the iteration loop
        unsigned int nNLIterations = 0;
        unsigned int nLSIterations = 0;
        float lambda = 1;
        bool converged, lsCheck;
        checkTolerance(R, tol, converged);
        unsigned int rank;        
        floatMatrix oldFloatOuts = floatOuts; //Copy the float outs
        intMatrix   oldIntOuts   = intOuts;   //Copy the int outs

        //Begin the iteration loop
        while ((!converged) && (nNLIterations<maxNLIterations)){

            //Perform the linear solve
            ddx = -vectorTools::solveLinearSystem( J, R, rank);
            dx += ddx;

            //Check the rank to make sure the linear system has a unique solution
            if (rank != R.size()){
                return new errorNode("newtonRaphson", "The jacobian matrix is singular");
            }

            //Compute the new residual
            error = residual(x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts);

            if (error){
                errorOut result = new errorNode("newtonRaphson", "Error in residual calculation in non-linear iteration"); 
                result->addNext(error);
                return result;
            }

            //Check the line search criteria
            checkLSCriteria(R, Rp, lsCheck, alpha);
            nLSIterations = 0;
            lambda = 1;

            //Enter line-search if required
            while ((!lsCheck) && (nLSIterations < maxLSIterations)){

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
                error = residual(x0 + dx, floatArgs, intArgs, R, J, floatOuts, intOuts);

                if (error){
                    errorOut result = new errorNode("newtonRaphson", "Error in line-search");
                    result->addNext(error);
                    return result;
                }

                //Perform the line-search check
                checkLSCriteria(R, Rp, lsCheck, alpha);

                //Increment the number of line-search iterations
                nLSIterations++;
            }

            if (!lsCheck){
                convergeFlag = false;
                return new errorNode("newtonRaphson", "The line-search failed to converge.");
            }
            else{
                Rp = R;
                oldFloatOuts = floatOuts;
                oldIntOuts   = intOuts;
            }

            //Check if the solution is converged
            checkTolerance(R, tol, converged);

            //Increment nNLIterations
            nNLIterations++;
        }

        //Check if the solution converged
        if (!converged){
            convergeFlag = false;
            return new errorNode("newtonRaphson", "The Newton-Raphson solver failed to converge.");
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
                            floatVector &x, bool &convergeFlag, floatMatrix &floatOuts, intMatrix &intOuts, 
                            const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola,
                            const floatType alpha, const unsigned int maxLSIterations, const unsigned int homotopySteps){
        /*!
         * Solve a non-linear equation using a homotopy Newton solver. This method 
         * can be successful in solving very stiff equations which other techniques 
         * struggle to capture. It effectively breaks the solve up into sub-steps 
         * of easier to solve equations which will eventually converge to the 
         * more difficult problem.
         * 
         * The residual function should have inputs of the form
         * :param const floatVector &x: A vector of the variable to be solved.
         * :param const floatMatrix &floatArgs: Additional floating point arguments to residual
         * :param const intMatrix &intArgs: Additional integer arguments to the residual
         * :param floatVector &residual: The residual vector
         * :param floatMatrix &jacobian: The jacobian matrix of the residual w.r.t. x
         * :param floatMatrix &floatOuts: Additional floating point values to return.
         * :param intMatrix &intOuts: Additional integer values to return.
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param floatVector &x: The converged value of the solver.
         * :param bool &convergeFlag: A flag which indicates whether the solver converged.
         * :param floatMatrix &floatOuts: Additional floating point values to return.
         * :param intMatrix &intOuts: Additional integer values to return.
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         * :param const unsigned int maxNLIterations: The maximum number of non-linear iterations.
         * :param const floatType tolr: The relative tolerance
         * :param const floatType tola: The absolute tolerance 
         * :param const floatType alpha: The line search criteria.
         * :param const unsigned int maxLSIterations: The maximum number of line-search iterations.
         * :param const unsigned int homotopySteps: The number of homotopy steps which will be taken.
         */

        //Initialize the homotopy solver
        floatType ds = 1./homotopySteps;
        floatType s  = 0;
        floatVector xh = x0;

        //Define the homotopy residual equation
        floatVector Rinit;
        floatMatrix J;

        //Compute the initial residual
        errorOut error = residual(x0, floatArgs, intArgs, Rinit, J, floatOuts, intOuts);

        if (error){
            errorOut result = new errorNode("homotopySolver", "error in initial residual calculation");
            result->addNext(error);
            return result;
        }

        //Define the homotopy residual
        stdFncNLFJ homotopyResidual;
        homotopyResidual = [&](const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                            floatVector &r, floatMatrix &J, floatMatrix &fO, intMatrix &iO){
 
            floatVector R;
            error = residual(x_, floatArgs_, intArgs_, R, J, fO, iO);
 
            if (error){
                errorOut result = new errorNode("homotopySolver::homotopyResidual", "error in residual calculation");
                result->addNext(error);
                return result;
            }
 
            r = R - (1 - s)*Rinit;
 
            return static_cast<errorOut>(NULL);
        };

        //Begin the homotopy loop
        for (unsigned int n=0; n<homotopySteps; n++){

            //Update s
            s += ds;

            error = newtonRaphson( homotopyResidual, xh, x, convergeFlag, floatOuts, intOuts, 
                                   floatArgs, intArgs, maxNLIterations, tolr, tola, 
                                   alpha, maxLSIterations);

            if (error){
                errorOut result = new errorNode("homotopySolver", "error in Newton Raphson solution");
                result->addNext(error);
                return result;
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
         * :param const floatVector &R: The residual vector.
         * :param const floatVector &tol: The tolerance.
         * :param bool result: The result
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
         * l2norm(R) < (1 - alpha)* l2norm(Rp)
         * 
         * :param const floatVector &R: The trial residual.
         * :param const floatVector &Rp: the previous acceptable residual
         * :param bool &result: The output value.
         * :param const floatType &alpha: The scaling factor on Rp
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
         * :param const floatVector &x: A vector of the variable to be solved.
         * :param const floatMatrix &floatArgs: Additional floating point arguments to function
         * :param const intMatrix &intArgs: Additional integer arguments to the function
         * :param floatVector &value: The output value vector
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param floatVector &grad: The finite difference gradient. 
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         * :param const floatType eps: The perturbation. delta[i] = eps*(x0[i]) + eps
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
         * :param const floatVector &x: A vector of the variable to be solved.
         * :param const floatMatrix &floatArgs: Additional floating point arguments to residual
         * :param const intMatrix &intArgs: Additional integer arguments to the residual
         * :param floatVector &residual: The residual vector
         * :param floatMatrix &jacobian: The jacobian matrix
         * :param floatMatrix &floatOuts: Additional returning floating point values.
         * :param int Matrix &intOuts: Additional return integer values.
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         * :param bool isGood: Whether the error in the jacobian is within tolerance.
         * :param const floatType eps: The perturbation. delta[i] = eps*(x0[i]) + eps
         * :param const floatType tolr: The relative tolerance
         * :param const floatType tola: The absolute tolerance
         * :param const bool suppressOutput: Suppress the output to the terminal
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
}
