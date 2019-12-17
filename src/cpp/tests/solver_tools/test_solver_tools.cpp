//Tests for constitutive_tools

#include<solver_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

typedef solverTools::errorOut errorOut;
typedef solverTools::errorNode errorNode;
typedef solverTools::floatType floatType;
typedef solverTools::floatVector floatVector;
typedef solverTools::floatMatrix floatMatrix;
typedef solverTools::intVector intVector;
typedef solverTools::intMatrix intMatrix;

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer)
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

errorOut nlFxn1(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                floatVector &residual, floatMatrix &jacobian){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * :param const floatVector &x: The variable vector
     * :param const floatMatrix &floatArgs: Floating point arguments to the function
     * :param const intMatrix &intArgs: Integer arguments to the function
     * :param floatVector &residual: The residual vector output.
     * :param floatMatrix &jacobian: The jacobian output.
     */

    if (x.size() != 2){
        return new errorNode("nlFnx1", "x must have a size of 2");
    }

    floatType x0 = -1;
    floatType y0 = 5.6;

    residual = {x[0] - x0, x[1] - y0};
    jacobian = {{1, 0}, {0, 1}};

    return NULL;
}

errorOut nlFxn1(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                floatVector &residual){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * :param const floatVector &x: The variable vector
     * :param const floatMatrix &floatArgs: Floating point arguments to the function
     * :param const intMatrix &intArgs: Integer arguments to the function
     * :param floatVector &residual: The residual vector output.
     */

    floatMatrix Jtmp;
    return nlFxn1(x, floatArgs, intArgs, residual, Jtmp);
}

errorOut nlFxn2(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                floatVector &residual, floatMatrix &jacobian){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * :param const floatVector &x: The variable vector
     * :param const floatMatrix &floatArgs: Floating point arguments to the function
     * :param const intMatrix &intArgs: Integer arguments to the function
     * :param floatVector &residual: The residual vector output.
     * :param floatMatrix &jacobian: The jacobian output.
     */

    if (x.size() != 3){
        return new errorNode("nlFnx1", "x must have a size of 2");
    }

    floatType x0 = -1;
    floatType y0 = 5.6;
    floatType z0 = 9.3;

    residual = {(x[0] - 1)*(x[0] - 7)*x[1], (x[1] - 1)*(x[0] - 3)*x[2], x[0]*x[1]*x[2]};
    jacobian = {{(x[0] - 7)*x[1] + (x[0] - 1)*x[1], (x[0] - 1)*(x[0] - 7), 0},
                {   (x[1] - 1)*x[2],    (x[0] - 3)*x[2], (x[1] - 1)*(x[0] - 3)},
                {   x[1]*x[2],    x[0]*x[2], x[0]*x[1]}};

    return NULL;
}

errorOut nlFxn2(const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                floatVector &residual){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * :param const floatVector &x: The variable vector
     * :param const floatMatrix &floatArgs: Floating point arguments to the function
     * :param const intMatrix &intArgs: Integer arguments to the function
     * :param floatVector &residual: The residual vector output.
     */

    floatMatrix Jtmp;
    return nlFxn2(x, floatArgs, intArgs, residual, Jtmp); 
}

int testCheckTolerance(std::ofstream &results){
    /*!
     * Test the tolerance checking function.
     * 
     * :param std::ofstream &results: The output file
     */

    floatVector R   = {  1,   2, 3.00000, -4.0};
    floatVector tol = {1.5, 2.1, 3.00001,  4.1};
    bool result;

    errorOut error = solverTools::checkTolerance(R, tol, result);

    if (error){
        error->print();
        results << "testCheckTolerance & False\n";
        return 1;
    }

    if (!result){
        results << "testCheckTolerance (test 1) & False\n";
        return 1;
    }

    tol[0] = .98;

    error = solverTools::checkTolerance(R, tol, result);

    if (error){
        error->print();
        results << "testCheckTolerance & False\n";
        return 1;
    }

    if (result){
        results << "testCheckTolerance (test 2) & False\n";
    }

    tol[0] = 1.5;
    tol[3] = 3.8;
    error = solverTools::checkTolerance(R, tol, result);

    if (error){
        error->print();
        results << "testCheckTolerance & False\n";
        return 1;
    }

    if (result){
        results << "testCheckTolerance (test 3) & False\n";
        return 1;
    }

    tol = {1.6, 2.5};

    error = solverTools::checkTolerance(R, tol, result);

    if (!error){
        results << "testCheckTolerance (test 4) &False\n";
        return 1;
    }

    results << "testCheckTolerance & True\n";
    return 0;
}

int testNewtonRaphson(std::ofstream &results){
    /*!
     * Tests of the Newton-Raphson solver
     * 
     * :param std::ofstream &results: The output file
     */

    //The first test
    floatVector x0 = {1.5, 6};
    floatVector x;

    errorOut error = solverTools::newtonRaphson(nlFxn1, x0, x, {}, {});

    if (error){
        error->print();
        results << "testNewtonRaphson & False\n";
        return 1;
    }

    floatVector Rtmp;
    floatMatrix Jtmp;
    error = nlFxn1(x, {}, {}, Rtmp, Jtmp);

    if (error){
        error->print();
        results << "testNewtonRaphson & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(Rtmp, {0, 0})){
        results << "testNewtonRaphson (test 1) & False\n";
        return 1;
    }

    //The second test
    x0 = {1, 1, 1};
    error = solverTools::newtonRaphson(nlFxn2, x0, x, {}, {});

    if (error){
        error->print();
        results << "testNewtonRaphson & False\n";
        return 1;
    }

    error = nlFxn2(x, {}, {}, Rtmp, Jtmp);

    if (error){
        error->print();
        results << "testNewtonRaphson & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(Rtmp, {0, 0, 0})){
        results << "testNewtonRaphson (test 1) & False\n";
        return 1;
    }

//    if (!vectorTools::fuzzyEquals(x, {-1, 5.6})){
//        results << "testNewtonRaphson (test 1) & False\n";
//        return 1;
//    }


    results << "testNewtonRaphson & True\n";
    return 0;
}

int testFiniteDifference(std::ofstream &results){
    /*!
     * Test the finite difference jacobian calculator.
     * 
     * :param std::ofstream &results: The output file
     */

    //The first test
    floatVector x0 = {1.5, 6};
    floatMatrix J;

    solverTools::finiteDifference(nlFxn1, x0, J, {}, {});

    floatVector Rtmp;
    floatMatrix result;
    nlFxn1(x0, {}, {}, Rtmp, result);

    if (!vectorTools::fuzzyEquals(J, result)){
        results << "testFiniteDifference (test 1) & False\n";
        return 1;
    }

    //The second test
    x0 = {1, 1, 1};
    solverTools::finiteDifference(nlFxn2, x0, J, {}, {});
    nlFxn2(x0, {}, {}, Rtmp, result);

    if (!vectorTools::fuzzyEquals(J, result)){
        results << "testFiniteDifference (test 2) & False\n";
        return 1;
    }

    results << "testFiniteDifference & True\n";
    return 0;
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    testCheckTolerance(results);
    testNewtonRaphson(results);
    testFiniteDifference(results);

    //Close the results file
    results.close();

    return 0;
}
