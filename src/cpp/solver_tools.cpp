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
        unsigned int niter = 0;
        

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
}
