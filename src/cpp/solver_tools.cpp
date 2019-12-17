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

    errorOut newtonRaphson( errorOut (*residual)(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs, 
                                                 floatVector &residual, floatMatrix &jacobian),
                            const floatVector &x0,
                            floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                            const unsigned int maxNLIterations, const floatType tolr, const floatType tola){
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
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param floatVector &x: The converged value of the solver. 
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         * :param const unsigned int maxNLIterations: The maximum number of non-linear iterations.
         * :param floatType tolr: The relative tolerance
         * :param floatType tola: The absolute tolerance 
         */

        //Compute the initial residual and jacobian
        floatVector dx = floatVector(x0.size(), 0);
        floatVector ddx;
        floatVector R;
        floatMatrix J;

        errorOut error = residual(x0 + dx, floatArgs, intArgs, R, J);

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

        //Initialize variables required for the iteration loop
        int nNLIterations = 0;
        bool converged;
        checkTolerance(R, tol, converged);
        unsigned int rank;        

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
            error = residual(x0 + dx, floatArgs, intArgs, R, J);

            if (error){
                return new errorNode("newtonRaphson", "Error in residual calculation in non-linear iteration"); 
            }

            //Check if the solution is converged
            checkTolerance(R, tol, converged);

            //Increment nNLIterations
            nNLIterations++;
        }

        //Check if the solution converged
        if (!converged){
            return new errorNode("newtonRaphson", "The Newton-Raphson solver failed to converge.");
        }
        else{
            //Update x
            x = x0 + dx;
            //Solver completed successfully
            return NULL;        
        }
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

    errorOut finiteDifference(errorOut (*residual)(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                                                 floatVector &residual),
                            const floatVector &x0,
                            floatMatrix &J, const floatMatrix &floatArgs, const intMatrix &intArgs, const floatType eps){
        /*!
         * Perform a forward finite difference gradient solve of the provided residual equation.
         * Note that for residual functions that are more (or less) complex than this you may need to 
         * wrap the function.
         * 
         * The residual function should have inputs of the form
         * :param const floatVector &x: A vector of the variable to be solved.
         * :param const floatMatrix &floatArgs: Additional floating point arguments to residual
         * :param const intMatrix &intArgs: Additional integer arguments to the residual
         * :param floatVector &residual: The residual vector
         * 
         * The main routine accepts the following parameters:
         * :param const floatVector &x0: The initial iterate of x.
         * :param floatVector &J: The finite difference jacobian. 
         * :param const floatMatrix &floatArgs: The additional floating-point arguments.
         * :param const intMatrix &intArgs: The additional integer arguments.
         */

        //Initialize the first residual and the gradient
        floatVector R0, Ri;

        //Compute the first residual
        errorOut error = residual(x0, floatArgs, intArgs, R0);
        if (error){
            errorOut result = new errorNode("finiteDifference", "Error in initial residual calculation");
            result->addNext(error);
            return result;
        }

        J = floatMatrix(R0.size(), floatVector(x0.size(), 0));

        for (unsigned int i=0; i<x0.size(); i++){
            //Set the step size
            floatVector delta = floatVector(x0.size(), 0);
            delta[i] = eps*fabs(x0[i]) + eps;

            //Compute the residual after perturbation
            error = residual(x0 + delta, floatArgs, intArgs, Ri);
            if (error){
                errorOut result = new errorNode("finiteDifference", "Error in residual calculation");
                result->addNext(error);
                return result;
            }

            //Set the terms of the gradient
            for (unsigned int j=0; j<Ri.size(); j++){
                J[j][i] = (Ri[j] - R0[j])/delta[i];
            }
        }
        return NULL;
    }
}
