//Tests for solver_tools

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

errorOut nlFxn1( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts,
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    if ( x.size( ) != 2 ){
        return new errorNode( "nlFnx1", "x must have a size of 2" );
    }

    floatType x0 = -1;
    floatType y0 = 5.6;

    residual = { x[0] - x0, x[1] - y0 };
    jacobian = { { 1, 0 }, { 0, 1 } };

    floatOuts = { { -1 }, { -1, -2, -3 }, { 4, 5, 6 } };
    intOuts = { { 1, 2, 8 } };

    return NULL;
}

errorOut nlFxn1( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual ){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     */

    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    return nlFxn1( x, floatArgs, intArgs, residual, Jtmp, fO, iO );
}

errorOut nlFxn2( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts, 
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    if ( x.size( ) != 3 ){
        return new errorNode( "nlFnx2", "x must have a size of 3" );
    }

    residual = { ( x[ 0 ] - 1 ) * ( x[ 0 ] - 7 ) * x[ 1 ], ( x[ 1 ] - 1 ) * ( x[ 0 ] - 3 ) * x[ 2 ], x[ 0 ] * x[ 1 ] * x[ 2 ] };
    jacobian = { { ( x[ 0 ] - 7 ) * x[ 1 ] + ( x[ 0 ] - 1 ) * x[ 1 ], ( x[ 0 ] - 1 ) * ( x[ 0 ] - 7 ), 0 },
                 {   ( x[ 1 ] - 1 ) * x[ 2 ],    ( x[ 0 ] - 3 ) * x[ 2 ], ( x[ 1 ] - 1 ) * ( x[ 0 ] - 3 ) },
                 {   x[ 1 ] * x[ 2 ],    x[ 0 ] * x[ 2 ], x[ 0 ] * x[ 1 ] } };

    floatOuts = { { -1 }, { -1, -2, -3 }, { 4, 5, 6 } };
    intOuts = { { 1, 2, 8 } };

    return NULL;
}

errorOut nlFxn2( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual ){
    /*!
     * A non-linear function for use in testing the solver.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     */

    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    return nlFxn2( x, floatArgs, intArgs, residual, Jtmp, fO, iO );
}

errorOut nlFxn3( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts, 
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver which will 
     * require the use of the line-search algorithm.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    residual = { std::exp( -x[ 0 ] ) - 1 };
    jacobian = { { -std::exp( -x[ 0 ] ) } };
    floatOuts = { { -1 }, { -1, -2, -3 }, { 4, 5, 6 } };
    intOuts = { { 1, 2, 8 } };
    return NULL;
}

errorOut nlFxn3( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual ){
    /*!
     * A non-linear function for use in testing the solver which will 
     * require the use of the line-search algorithm.
     * 
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     */
    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    return nlFxn3( x, floatArgs, intArgs, residual, Jtmp, fO, iO );
}

errorOut nlFxn4( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts,
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver which will
     * require the use of the line-search algorithm.
     *
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    residual = { std::tanh( x[ 0 ] ) };
    jacobian = { { ( std::cosh( x[ 0 ] ) * std::cosh( x[ 0 ] ) - std::sinh( x[ 0 ] ) * std::sinh( x[ 0 ] ) ) / ( std::cosh( x[ 0 ] ) * std::cosh( x[ 0 ] ) ) } };
    return NULL;
}

errorOut nlFxn5( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts,
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver which will require
     * the use of the bounded homotopy solver
     *
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    //floatArgs answers
    floatVector answer1 = { .1, .2, .3, .4 };
    floatVector answer2 = { -0.01, -0.02 };

    //IntArgs answers
    intVector answer3 = { -1, -2, -3 };
    intVector answer4 = { 5, 4, 3, 2 };
    intVector answer5 = { 8, 9, 9 };

    //floatOuts answers
    floatVector answer6 = { 0, 1, 2 };
    floatVector answer7 = { 7, -6 };
    floatVector answer8 = { .24, .25 };

    //intOuts answers
    intVector   answer9  = { 1, 2, 3 };
    intVector   answer10 = { -5, 6, 7, 8 };

    //x tests
    if ( x.size() != 1 ){
        return new errorNode( "nlFxn5", "The x vector should have a size of 1" );
    }

    //floatArgs tests
    if ( floatArgs.size() != 2 ){
        return new errorNode( "nlFxn5", "The floatArgs matrix should have two values" );
    }

    if ( !vectorTools::fuzzyEquals( floatArgs[ 0 ], answer1) ){
        return new errorNode( "nlFxn5", "The first value of floatArgs should be { 0.1, 0.2, 0.3, 0.4 }" );
    }

    if ( !vectorTools::fuzzyEquals( floatArgs[ 1 ], answer2) ){
        return new errorNode( "nlFxn5", "The second value of floatArgs should be { -0.01, -0.02 }" );
    }

    //intArgs tests
    if ( intArgs.size() != 3 ){
        return new errorNode( "nlFxn5", "The intArgs matrix should have three values" );
    }

    if ( !vectorTools::fuzzyEquals( intArgs[ 0 ], answer3 ) ){
        return new errorNode( "nlFxn5", "The first value of intargs should be { -1, -2, -3 }" );
    }

    if ( !vectorTools::fuzzyEquals( intArgs[ 1 ], answer4 ) ){
        return new errorNode( "nlFxn5", "The second value of intargs should be { 5, 4, 3, 2 }" );
    }

    if ( !vectorTools::fuzzyEquals( intArgs[ 2 ], answer5 ) ){
        return new errorNode( "nlFxn5", "The third value of intargs should be { 8, 9, 9 }" );
    }

    //floatOuts tests
    if ( floatOuts.size() != 3 ){
        return new errorNode( "nlFxn5", "The floatOuts matrix should have three values" );
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 0 ], answer6 ) ){
        return new errorNode( "nlFxn5", "The first values in the floatOuts should be { 0, 1, 2 }" );
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 1 ], answer7 ) ){
        return new errorNode( "nlFxn5", "The second values in the floatOuts should be { 7, -6 }" );
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 2 ], answer8 ) ){
        return new errorNode( "nlFxn5", "The third values in the floatOuts should be { 0.24, 0.25 }" );
    }

    //intOuts tests
    if ( intOuts.size() != 2 ){
        return new errorNode( "nlFxn5", "The intOuts matrix must have a size of 2" );
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 0 ], answer9 ) ){
        return new errorNode( "nlFxn5", "The first values in the intOuts should be { 1, 2, 3 }" );
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 1 ], answer10 ) ){
        return new errorNode( "nlFxn5", "The second values in the intOuts should be { -5, 6, 7, 8 }" );
    }

    residual = { ( x[ 0 ] - 1. ) * ( x[ 0 ] + 1 ) * ( x[ 0 ] - 0.25 ) * ( x[ 0 ] + 0.1 ) };

    jacobian = { {  ( x[ 0 ] + 1. ) * ( x[ 0 ] - 0.25 ) * ( x[ 0 ] + 0.1  )
                  + ( x[ 0 ] - 1. ) * ( x[ 0 ] - 0.25 ) * ( x[ 0 ] + 0.1  )
                  + ( x[ 0 ] - 1. ) * ( x[ 0 ] + 1.   ) * ( x[ 0 ] + 0.1  )
                  + ( x[ 0 ] - 1. ) * ( x[ 0 ] + 1.   ) * ( x[ 0 ] - 0.25 ) } };

    floatOuts = { floatOuts[ 0 ] + 0.1, floatOuts[ 1 ], floatOuts[ 0 ] };
    intOuts = { intOuts[ 0 ] - 2, intOuts[ 0 ], intOuts[ 1 ] };

    return NULL;
}

errorOut nlFxn6( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts,
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver which will require
     * the use of the bounded homotopy solver
     *
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    floatType x1 = x[ 0 ];
    floatType x2 = x[ 1 ];
    floatType x3 = x[ 2 ];

    residual.resize( 3 );

    residual[ 0 ] = ( x1 - 1 )*( x1 + 1 )*( x1 - 0.25 )*( x1 + 0.1 );
    residual[ 1 ] = ( x2 - 1 ) * ( x2 - 1 );
    residual[ 2 ] = ( x1 + 5 ) * ( x3 + 1 );

    floatType dr1dx1 = ( x1 - 1 ) * ( x1 - 0.25 ) * ( x1 + 0.1 )
                     + ( x1 - 1 ) * ( x1 - 0.25 ) * ( x1 + 1 )
                     + ( x1 - 1 ) * ( x1 + 0.1 ) * ( x1 + 1 )
                     + ( x1 - 0.25 ) * ( x1 + 0.1 ) * ( x1 + 1 );

    floatType dr1dx2 = 0.;
    floatType dr1dx3 = 0.;

    floatType dr2dx1 = 0.;
    floatType dr2dx2 = 2 * ( x2 - 1 );
    floatType dr2dx3 = 0.;

    floatType dr3dx1 = x3 + 1;
    floatType dr3dx2 = 0.;
    floatType dr3dx3 = x1 + 5;

    jacobian = { { dr1dx1, dr1dx2, dr1dx3 },
                 { dr2dx1, dr2dx2, dr2dx3 },
                 { dr3dx1, dr3dx2, dr3dx3 } };

    return NULL;
}

errorOut nlFxn7( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual, floatMatrix &jacobian, floatMatrix &floatOuts,
                 intMatrix &intOuts ){
    /*!
     * A non-linear function for use in testing the solver which will require
     * the use of the bounded homotopy solver
     *
     * /param &x: The variable vector
     * /param &floatArgs: Floating point arguments to the function
     * /param &intArgs: Integer arguments to the function
     * /param &residual: The residual vector output.
     * /param &jacobian: The jacobian output.
     * /param &floatOuts: Additional floating point outputs.
     * /param &intOuts: Additional integer outputs.
     */

    if ( x.size() != 1 ){
        return new errorNode( "nlFxn7", "The x vector must have a size of 1" );
    }

    residual = { std::log( x[ 0 ] ) };

    jacobian = { { 1. / x[ 0 ] } };

    return NULL;
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

int testNewtonRaphson( std::ofstream &results ){
    /*!
     * Tests of the Newton-Raphson solver
     * 
     * /param &results: The output file
     */

    //The first test
    floatVector x0 = { 1.5, 6 };
    floatVector x;
    bool converged, fatalError;

    solverTools::stdFncNLFJ func;
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn1 );

    floatMatrix floatOut;
    intMatrix intOut;    
    errorOut error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { } );

    if ( error ){
        error->print( );
        results << "testNewtonRaphson nlFxn1 & False\n";
        return 1;
    }

    floatVector Rtmp;
    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    error = nlFxn1( x, { }, { }, Rtmp, Jtmp, fO, iO );

    if ( error ){
        error->print( );
        results << "testNewtonRaphson nlFxn1 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Rtmp, { 0, 0 } ) ){
        results << "testNewtonRaphson (test 1) & False\n";
        return 1;
    }

    //The second test
    x0 = { 1, 1, 1 };
    floatOut.clear( );
    intOut.clear( );
    fO.clear( );
    iO.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn2 );
    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { } );

    if ( error ){
        error->print( );
        results << "testNewtonRaphson nlFxn2 & False\n";
        return 1;
    }

    error = nlFxn2( x, { }, { }, Rtmp, Jtmp, fO, iO );

    if ( error ){
        error->print( );
        results << "testNewtonRaphson nlFxn2 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Rtmp, { 0, 0, 0 } ) ){
        results << "testNewtonRaphson (test 2) & False\n";
        return 1;
    }

    //The third test
    x0 = { 3 };
    floatOut.clear( );
    intOut.clear( );
    fO.clear( );
    iO.clear( );
    
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn3 );
    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { } );

    if ( error ){
        error->print( );
        results << "testNewtonRaphson nlFxn3 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x, { 0 } ) ){
        results << "testNewtonRaphson (test 3) & false\n";
        return 1;
    }

    //The fourth test. Tests the bounded Newton method
    x0 = { 10. };
    floatOut.clear();
    intOut.clear();
    fO.clear();
    iO.clear();

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn7 );

    solverTools::solverType linearSolver;
    floatMatrix J;

    intVector boundVariableIndices = { 0 };
    intVector boundSigns = { 0 };
    floatVector boundValues = { 1e-9 };
    floatMatrix Jexp = { { 1. } };

    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, {}, {}, linearSolver, J );

    if ( error ){
        error->print();
        results << "testNewtonRaphson nlFxn7 & False\n";
        return 1;
    }
    if ( !vectorTools::fuzzyEquals( x, { 1. } ) ){
        results << "testNewtonRaphson (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( J, Jexp ) ){
        results << "testNewtonRaphson (test 5) & False\n";
        return 1;
    }

    //The fifth test: This test makes sure that when a Newton-Raphson iteration fails it
    //correctly returns a failure to converge.
    x0 = { -5. };
    floatOut.clear( );
    intOut.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn4 );

    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { }, 5 );

    if ( converged ){
        results << "testNewtonRaphson (test 6) & False\n";
        return 1;
    }

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

    solverTools::stdFncNLF func;
    func = static_cast<solverTools::NonLinearFunction>(nlFxn1);
    solverTools::finiteDifference(func, x0, J, {}, {});

    floatVector Rtmp;
    floatMatrix result;
    floatMatrix floatOuts;
    intMatrix intOuts;
    nlFxn1(x0, {}, {}, Rtmp, result, floatOuts, intOuts);

    if (!vectorTools::fuzzyEquals(J, result)){
        results << "testFiniteDifference (test 1) & False\n";
        return 1;
    }

    //The second test
    x0 = {1, 1, 1};
    floatOuts.clear();
    intOuts.clear();
    func = static_cast<solverTools::NonLinearFunction>(nlFxn2);
    solverTools::finiteDifference(func, x0, J, {}, {});
    nlFxn2(x0, {}, {}, Rtmp, result, floatOuts, intOuts);

    if (!vectorTools::fuzzyEquals(J, result)){
        results << "testFiniteDifference (test 2) & False\n";
        return 1;
    }

    results << "testFiniteDifference & True\n";
    return 0;
}

int testCheckJacobian(std::ofstream &results){
    /*!
     * Test the jacobian checking utility.
     * 
     * :param std::ofstream &results: The output file
     */

    //The first test

    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>(nlFxn1);
    bool isGood = false;
    floatType eps = 1e-6;
    floatType tolr = 1e-6;
    floatType tola = 1e-6;
    bool suppressOutput = true;

    floatVector x0 = {0, 0};
    errorOut error = solverTools::checkJacobian(func, x0, {}, {}, isGood, eps, tolr, tola, suppressOutput);

    if (error){
        error->print();
        results << "testCheckJacobian & False\n";
        return 1;
    }

    if (!isGood){
        results << "testCheckJacobian (test 1) & False\n";
        return 1;
    }

    //The second test
    solverTools::stdFncNLFJ badfunc;
    badfunc = [&](const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                            floatVector &r, floatMatrix &j, floatMatrix &fO, intMatrix &iO){
        errorOut e = func(x_, floatArgs_, intArgs_, r, j, fO, iO);
        j[0][1] = 0.1;
        return e;
    }; 

    error = solverTools::checkJacobian(badfunc, x0, {}, {}, isGood, eps, tolr, tola, suppressOutput);

    if (error){
        error->print();
        results << "testCheckJacobian & False\n";
        return 1;
    }

    if (isGood){
        results << "testCheckJacobian (test 2) & False\n";
        return 1;
    }

    results << "testCheckJacobian & True\n";
    return 0;
}

int testCheckLSCriteria(std::ofstream &results){
    /*!
     * Test the line search criteria
     * 
     * :param std::ofstream &results: The output file
     */

    floatVector R  = {1, 2, 3, 4, 5, 6};
    floatVector Rp = {2, 3, 4, 5, 6, 7};
    bool result;

    solverTools::checkLSCriteria(R, Rp, result);

    if (!result){
        results << "testCheckLSCriteria (test 1) & False\n";
        return 1;
    }

    R[0] = 100;

    solverTools::checkLSCriteria(R, Rp, result);

    if (result){
        results << "testCheckLSCriteria (test 2) & False\n";
        return 1;
    }

    results << "testCheckLSCriteria  & True\n";
    return 0;
}

int testHomotopySolver(std::ofstream &results){
    /*!
     * Test the Homotopy solver.
     * 
     * :param std::ofstream &results: The output file
     */

    //The first test
    floatVector x0 = { 1.5, 6 };
    floatVector x;
    bool converged, fatalErrorFlag;
    floatMatrix floatOuts;
    intMatrix intOuts;

    solverTools::stdFncNLFJ func;
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn1 );
   
    errorOut error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { } );

    if ( error ){
        error->print( );
        results << "testHomotopySolver & False\n";
        return 1;
    }

    floatVector Rtmp;
    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    error = nlFxn1( x, { }, { }, Rtmp, Jtmp, fO, iO );

    if ( error ){
        error->print( );
        results << "testHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Rtmp, { 0, 0 } ) ){
        std::cout << "Rtmp: "; vectorTools::print( Rtmp );
        results << "testHomotopySolver (test 1) & False\n";
        return 1;
    }

    //The second test
    x0 = {1, 1, 1};
    floatOuts.clear();
    intOuts.clear();
    fO.clear();
    iO.clear();

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn2 );
    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { } );

    if ( error ){
        error->print( );
        results << "testHomotopySolver & False\n";
        return 1;
    }

    error = nlFxn2( x, { }, { }, Rtmp, Jtmp, fO, iO );

    if ( error ){
        error->print( );
        results << "testHomotpySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( Rtmp, { 0, 0, 0 } ) ){
        results << "testHomotopySolver (test 2) & False\n";
        return 1;
    }

    //The third test
    x0 = { 3 };
    floatOuts.clear( );
    intOuts.clear( );
    fO.clear( );
    iO.clear( );
    
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn3 );
    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { },
                                         20, 1e-9, 1e-9, 1e-4, 5, 0.2, 0.01 );

    if ( error ){
        error->print( );
        results << "testHomotopySolver & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals( x, { 0 } ) ){
        results << "testHomotopySolver (test 3) & false\n";
        return 1;
    }

    //The fourth test
    x0 = { 3 };
    floatOuts.clear( );
    intOuts.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn4 );

    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { },
                                         20, 1e-9, 1e-9, 1e-4, 4, 1.0, 0.1 );

    if ( error ){
        error->print();
        results << "testHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x, { 0 } ) ){
        results << "testHomotopySolver (test 4) & False\n";
        return 1;
    }

    results << "testHomotopySolver & True\n";
    return 0;
    
}

int test_applyBoundaryLimitation( std::ofstream &results ){
    /*!
     * Test of the application of the boundary conditions.
     *
     * \param &results: The output file
     */

    floatVector x0 = { 1.0,  2.0, 3.0, -1.0, -2.0, -3.0 };
    floatVector dxDefault = { 0.1, -0.5, 1.0,  0.1,  2.1, -0.5 };

    intVector variableIndices = {    0,    3, 4   };
    intVector barrierSigns    = {    1,    0, 1   };
    floatVector barrierValues = { 1.05, -1.0, 0.0 };

    floatVector dx = dxDefault;

    floatVector xAnswer1 = { 1.05, 1.75, 3.5, -0.95, -0.95, -3.25 };
    floatVector xAnswer2 = { 1.033333, 1.833333, 3.333333, -1, -1.3, -3.1666667 };
    floatVector xAnswer3 = { 1.05, 1.5, 4, -0.8, 0.0, -3.5 };

    errorOut error = solverTools::applyBoundaryLimitation( x0, variableIndices, barrierSigns, barrierValues, dx );

    if ( error ){
        error->print();
        results << "test_applyBoundaryLimitation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x0 + dx, xAnswer1 ) ){
        results << "test_applyBoundaryLimitation (test 1) & False\n";
        return 1;
    }

    dx = dxDefault;
    x0[ 3 ] = -0.9;
    dx[ 3 ] = -0.3;

    error = solverTools::applyBoundaryLimitation( x0, variableIndices, barrierSigns, barrierValues, dx );

    if ( error ){
        error->print( );
        results << "test_applyBoundaryLimitation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x0 + dx, xAnswer2 ) ){
        results << "test_applyBoundaryLimitation (test 2) & False\n";
        return 1;
    }

    dx = dxDefault;
    error = solverTools::applyBoundaryLimitation( x0, variableIndices, barrierSigns, barrierValues, dx, 1e-9, 1e-9, true );

    if ( error ){
        error->print( );
        results << "test_applyBoundaryLimitation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x0 + dx, xAnswer3 ) ){
        results << "test_applyBoundaryLimitation (test 3) & False\n";
        return 1;
    }

    results << "test_applyBoundaryLimitation & True\n";
    return 0;
}

int main( ){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open( "results.tex" );

    //Run the tests
    testCheckTolerance( results );
    testNewtonRaphson( results );
    testFiniteDifference( results );
    testCheckJacobian( results );
    testCheckLSCriteria( results );
    testHomotopySolver( results );
    test_applyBoundaryLimitation( results );

    //Close the results file
    results.close( );

    return 0;
}
