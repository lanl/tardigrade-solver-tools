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

    //The fifth test ( hard bounds )
    x0 = { 10 };
    floatOuts.clear();
    intOuts.clear();

    solverTools::solverType linearSolver;
    floatMatrix J, Jexp;

    Jexp = { { 1 } };

    intVector variableIndices = { 0 };
    intVector barrierSigns = { 0 };
    floatVector barrierValues = { 1e-9 };

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn7 );

    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, {}, {},
                                         linearSolver, J, variableIndices, barrierSigns, barrierValues,
                                         false, 20, 1e-9, 1e-9, 1e-4, 4, 1.0, 0.1 );

    if ( error ){
        error->print();
        results << "testHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x, { 1 } ) ){
        results << "testHomotopySolver (test 5) & false\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( J, Jexp ) ){
        results << "testHomotopySolver (test 6) & False\n";
        return 1;
    }

    results << "testHomotopySolver & True\n";
    return 0;
    
}

int test_aFxn( std::ofstream &results ){
    /*!
     * Test the computation of the a parameter in the Barrier function.
     *
     * :param std::ofstream &results: The output file.
     */

    floatType pseudoT = .72;
    floatType logAfxn = 5.2;

    floatType answer  = 42.26671935907283;

    floatType result;

    errorOut error = solverTools::aFxn( pseudoT, logAfxn, result );

    if ( error ){
        error->print();
        results << "test_aFxn & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_aFxn (test 1) & False\n";
        return 1;
    }

    floatType dadT;

    error = solverTools::aFxn( pseudoT, logAfxn, result, dadT );

    if ( error ){
        error->print();
        results << "test_aFxn & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_aFxn (test 2) & False\n";
        return 1;
    }

    floatType eps = 1e-6;
    floatType delta = eps * pseudoT + eps;

    floatType aP, aM;

    error = solverTools::aFxn( pseudoT + delta, logAfxn, aP );

    if ( error ){
        error->print();
        results << "test_aFxn & False\n";
        return 1;
    }

    error = solverTools::aFxn( pseudoT - delta, logAfxn, aM );

    if ( error ){
        error->print();
        results << "test_aFxn & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( aP - aM ) / ( 2 * delta ), dadT ) ){
        results << "test_aFxn (test 3) & False\n";
        return 1;
    }

    results << "test_aFxn & True\n";
    return 0;
}

int test_computeBarrierFunction( std::ofstream &results ){
    /*!
     * Test the computation of the boundary function
     *
     * :param std::ofstream &results: The output file.
     */

    floatType x        = 0.4;
    floatType pseudoT  = 0.25;
    floatType logAmax  = 5;
    floatType b        = 0.14;

    floatType negativeSignAnswer = -0.5964638357684787;
    floatType positiveSignAnswer = 1.4780926435784547;

    floatType result;

    errorOut error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, false, result );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, negativeSignAnswer ) ){
        results << "test_computeBarrierFunction (test 1) & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, true, result );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, positiveSignAnswer ) ){
        results << "test_computeBarrierFunction (test 2) & False\n";
        return 1;
    }

    //Test the Jacobians
    floatType dbdx, dbdt;

    //Test the Jacobians when the sign is negative
    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, false, result, dbdx, dbdt );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, negativeSignAnswer ) ){
        results << "test_computeBarrierFunction (test 3) & False\n";
        return 1;
    }

    floatType eps = 1e-6;

    floatType dx = eps * fabs( x ) + eps;
    floatType dt = eps * fabs( pseudoT ) + eps;

    floatType bP, bM;

    error = solverTools::computeBarrierFunction( x + dx, pseudoT, logAmax, b, false, bP );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x - dx, pseudoT, logAmax, b, false, bM );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dx ), dbdx ) ){
        results << "test_computeBarrierFunction (test 4) & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x, pseudoT + dt, logAmax, b, false, bP );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x, pseudoT - dt, logAmax, b, false, bM );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dt ), dbdt ) ){
        results << "test_computeBarrierFunction (test 5) & False\n";
        return 1;
    }

    //Test the Jacobians when the sign is positive
    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, true, result, dbdx, dbdt );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, positiveSignAnswer ) ){
        results << "test_computeBarrierFunction (test 6) & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x + dx, pseudoT, logAmax, b, true, bP );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x - dx, pseudoT, logAmax, b, true, bM );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dx ), dbdx ) ){
        results << "test_computeBarrierFunction (test 7) & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x, pseudoT + dt, logAmax, b, true, bP );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    error = solverTools::computeBarrierFunction( x, pseudoT - dt, logAmax, b, true, bM );

    if ( error ){
        error->print();
        results << "test_computeBarrierFunction & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dt ), dbdt ) ){
        results << "test_computeBarrierFunction (test 8) & False\n";
        return 1;
    }

    results << "test_computeBarrierFunction & True\n";
    return 0;
}

int test_computeBarrierHomotopyResidual( std::ofstream &results ){
    /*!
     * Test the computation of the barrier homotopy residual
     *
     * :param std::ofstream &results: The output file
     */


    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>(nlFxn5);

    floatVector x0 = { 0. };
    floatMatrix floatArgsDefault =
        {
            { 0.28 }, //The pseudo-time
            { 0.1 },  //The barrier value
            { 10. },  //The logAMax values
            { .1, .2, .3, .4 }, //Function Parameters
            { -0.01, -0.02 }
        };

    intMatrix intArgsDefault =
        {
            { 0 }, //Variable indices
            { 0 }, //Residual indices
            { 0 }, //Barrier signs
            { -1, -2, -3 }, //Function parameters
            { 5, 4, 3, 2 },
            { 8, 9, 9 }
        };

    floatVector residualResult;
    floatMatrix jacobian;

    floatMatrix floatOutsDefault =
        {
            { 0, 1, 2 },
            { 7, -6 },
            { .24, .25 }
        };

    intMatrix intOutsDefault =
        {
            { 1, 2, 3 },
            { -5, 6, 7, 8 },
        };

    floatMatrix floatArgs = floatArgsDefault;
    intMatrix   intArgs   = intArgsDefault;
    floatMatrix floatOuts = floatOutsDefault;
    intMatrix   intOuts   = intOutsDefault;

    floatVector residualAnswer = { 0.2775586103363596 };

#ifdef DEBUG_MODE
    debugMap DEBUG;
#endif

    errorOut error = solverTools::computeBarrierHomotopyResidual( func, x0, floatArgs, intArgs, residualResult, jacobian,
                                                                  floatOuts, intOuts
#ifdef DEBUG_MODE
                                                                  , DEBUG
#endif
                                                                );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( residualAnswer, residualResult ) ){
        results << "test_computeBarrierHomotopyResidual (test 1) & False\n";
        return 1;
    }

    //Check that the non-homotopy outputs are as expected.
    if ( floatOuts.size() != 4 ){
        results << "test_computeBarrierHomotopyResidual (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 1 ], floatOutsDefault[ 0 ] + 0.1 ) ){
        results << "test_computeBarrierHomotopyResidual (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 2 ], floatOutsDefault[ 1 ] ) ){
        results << "test_computeBarrierHomotopyResidual (test 4) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 3 ], floatOutsDefault[ 0 ] ) ){
        results << "test_computeBarrierHomotopyResidual (test 5) & False\n";
        return 1;
    }

    if ( intOuts.size() != 3 ){
        results << "test_computeBarrierHomotopyResidual (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 0 ], intOutsDefault[ 0 ] - 2 ) ){
        results << "test_computeBarrierHomotopyResidual (test 7) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 1 ], intOutsDefault[ 0 ] ) ){
        results << "test_computeBarrierHomotopyResidual (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 2 ], intOutsDefault[ 1 ] ) ){
        results << "test_computeBarrierHomotopyResidual (test 9) & False\n";
        return 1;
    }

    //Test the Jacobians
    floatType eps = 1e-6;

    //Test drdx
    floatVector dx = eps * x0 + eps;

    floatVector rP, rM;
    floatMatrix JP, JM;

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0 + dx, floatArgs, intArgs, rP, JP,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual & False\n";
        return 1;
    }

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0 - dx, floatArgs, intArgs, rM, JM,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual & False\n";
        return 1;
    }

    floatVector gradCol = ( rP - rM ) / ( 2 * dx[ 0 ] );

    if ( !vectorTools::fuzzyEquals( gradCol, jacobian[ 0 ] ) ){
        results << "test_computeBarrierHomotopyResidual (test 10) & False\n";
        return 1;
    }

    //test drdt
    eps = 1e-7;
    floatType dt = eps * floatArgsDefault[ 0 ][ 0 ] + eps;

    floatArgs = floatArgsDefault;
    floatArgs[ 0 ][ 0 ] += dt;

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0, floatArgs, intArgs, rP, JP,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual & False\n";
        return 1;
    }

    floatArgs = floatArgsDefault;
    floatArgs[ 0 ][ 0 ] -= dt;

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0, floatArgs, intArgs, rM, JM,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual & False\n";
        return 1;
    }

    gradCol = ( rP - rM ) / ( 2 * dt );

    if ( !vectorTools::fuzzyEquals( gradCol, floatOuts[ 0 ], 1e-5 ) ){
        results << "test_computeBarrierHomotopyResidulal (test 11) & False\n";
        return 1;
    }

    results << "test_computeBarrierHomotopyResidual & True\n";
    return 0;
}

int test_computeBarrierHomotopyResidual2( std::ofstream &results ){
    /*!
     * Test for the computation of the barrier homotopy residual.
     *
     * :param std::ofstream &results: The output file.
     */


    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>(nlFxn6);

    floatVector          x = { 0.15, 0.1, -1.2 };
    floatType   pseudoTime = 0.24;

    floatVector logAMaxVals = { 10, 6 };
    floatVector bvals       = { 0.1, -1.1 };

    intVector   variableIndices = { 0, 2 };
    intVector   residualIndices = { 2, 1 };
    intVector   signs           = { 0, 1 };

    floatMatrix floatArgs = { { pseudoTime }, bvals, logAMaxVals };
    intMatrix   intArgs   = { variableIndices, residualIndices, signs };

    floatMatrix floatOuts = { {} };
    intMatrix   intOuts   = { {} };

    floatVector residualAnswer = { 0.0244375 ,  0.53651154, -0.97499937 };

    floatVector residualResult;
    floatMatrix jacobian;

    errorOut error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, residualResult, jacobian,
                                                                  floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual2 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( residualResult, residualAnswer ) ){
        results << "test_computeBarrierHomotopyResidual2 (test 1) & False\n";
        return 1;
    }

    //Tests of the Jacobians

    //Test the Jacobian w.r.t. x
    floatType eps = 1e-7;
    for ( unsigned int i = 0; i < x.size(); i++ ){
        floatVector delta( x.size(), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        floatVector rP, rM;

        error = solverTools::computeBarrierHomotopyResidual( func, x + delta, floatArgs, intArgs, rP, jacobian,
                                                             floatOuts, intOuts );

        if ( error ){
            error->print();
            results << "test_computeBarrierHomotopyResidual2 & False\n";
            return 1;
        }

        error = solverTools::computeBarrierHomotopyResidual( func, x - delta, floatArgs, intArgs, rM, jacobian,
                                                             floatOuts, intOuts );

        if ( error ){
            error->print();
            results << "test_computeBarrierHomotopyResidual2 & False\n";
            return 1;
        }

        floatVector gradCol = ( rP - rM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) ){
                results << "test_computeBarrierHomotopyResult2 (test 2) & False\n";
                return 1;
            }
        }
    }

    //Test the Jacobian w.r.t. t
    floatType dt = eps * fabs( pseudoTime ) + eps;

    floatVector rP, rM;

    floatArgs[ 0 ][ 0 ] = pseudoTime + dt;

    error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, rP, jacobian,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual2 & False\n";
        return 1;
    }

    floatArgs[ 0 ][ 0 ] = pseudoTime - dt;

    error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, rM, jacobian,
                                                         floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_computeBarrierHomotopyResidual2 & False\n";
        return 1;
    }

    floatVector gradCol = ( rP - rM ) / ( 2 * dt );

    if ( !vectorTools::fuzzyEquals( gradCol, floatOuts[ 0 ] ) ){
        results << "test_computeBarrierHomotopyResidual2 (test 3) & False\n";
        return 1;
    }

    results << "test_computeBarrierHomotopyResidual2 & True\n";
    return 0;
}

int test_barrierHomotopySolver( std::ostream &results ){
    /*!
     * Test the barrier Homotopy solver. This solver enables the addition of
     * bounds to a non-linear solve which can help prevent solutions from being
     * pulled into undesirable domains without having to resort to computing the
     * Hessian of the residual function as would be required for optimization
     * based techniques.
     *
     * :param std::ofstream &results: The output file.
     */

    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>(nlFxn5);

    floatVector barrierValues = { 0.1 };
    floatVector logAMaxValues = { 10. };
    floatMatrix floatArgsDefault =
        {
            { .1, .2, .3, .4 },
            { -0.01, -0.02 }
        };

    intVector variableIndices = { 0 };
    intVector residualIndices = { 0 };
    intVector barrierSigns    = { 0 };

    intMatrix intArgsDefault =
        {
            { -1, -2, -3 },
            { 5, 4, 3, 2 },
            { 8, 9, 9 }
        };

    floatVector residualResult;
    floatMatrix jacobian;

    floatMatrix floatOutsDefault =
        {
            { 0, 1, 2 },
            { 7, -6 },
            { .24, .25 }
        };

    intMatrix intOutsDefault =
        {
            { 1, 2, 3 },
            { -5, 6, 7, 8 },
        };

    floatMatrix floatArgs = floatArgsDefault;
    intMatrix   intArgs   = intArgsDefault;
    floatMatrix floatOuts = floatOutsDefault;
    intMatrix   intOuts   = intOutsDefault;

    floatType dt = 0.1;
    floatVector x0 = { 0. };
    bool implicitRefine = false;

    bool convergeFlag, fatalErrorFlag;

    floatVector result;
    floatVector answer = { 0.25 };

    errorOut error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                         barrierValues, logAMaxValues, floatArgs, intArgs,
                                                         implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                         floatOuts, intOuts,
                                                         20, 1e-9, 1e-9, 1e-4, 5, true );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_barrierHomotopySolver (test 1) & False\n";
        return 1;
    }

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    solverTools::solverType linearSolver;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                         barrierValues, logAMaxValues, floatArgs, intArgs,
                                                         implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                         floatOuts, intOuts, linearSolver, jacobian,
                                                         20, 1e-9, 1e-9, 1e-4, 5, true );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_barrierHomotopySolver (test 2) & False\n";
        return 1;
    }

    floatMatrix jacobianResult;

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = func( result, floatArgs, intArgs, residualResult, jacobianResult, floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( jacobian, jacobianResult ) ){
        results << "test_barrierHomotopySolver (test 3) & False\n";
        return 1;
    }

    x0 = { 0. };
    implicitRefine = true;

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                barrierValues, logAMaxValues, floatArgs, intArgs,
                                                implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts,
                                                20, 1e-9, 1e-9, 1e-4, 5, true );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_barrierHomotopySolver (test 4) & False\n";
        return 1;
    }

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                barrierValues, logAMaxValues, floatArgs, intArgs,
                                                implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts, linearSolver, jacobian,
                                                20, 1e-9, 1e-9, 1e-4, 5, true );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_barrierHomotopySolver (test 5) & False\n";
        return 1;
    }

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = func( result, floatArgs, intArgs, residualResult, jacobianResult, floatOuts, intOuts );

    if ( error ){
        error->print();
        results << "test_barrierHomotopySolver & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( jacobian, jacobianResult ) ){
        results << "test_barrierHomotopySolver (test 6) & False\n";
        return 1;
    }

    results << "test_barrierHomotopySolver & True\n";
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
        results << "test_applyBoundaryLimitation (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x0 + dx, xAnswer3 ) ){
        results << "test_applyBoundaryLimitation (test 3) & False\n";
        return 1;
    }

    results << "test_applyBoundaryLimitation & True\n";
    return 0;
}

int test_BFGS( std::ofstream &results ){
    /*!
     * Test of the BFGS optimization algorithm.
     *
     * :param std::ofstream &results: The output file
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>(lagrangian1);

    solverTools::floatVector x0 = { 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatVector xAnswer = { -1 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3}, {-0.4, -0.5, -0.6 } };
    intMatrix intOutsAnswer = { { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::BFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                        floatOuts, intOuts, floatArgs, intArgs );

    if ( error ){
        error->print();
        results << "test_BFGS & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( x, xAnswer ) ){
        results << "test_BFGS (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) ){
        results << "test_BFGS (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) ){
        results << "test_BFGS (test 3) & False\n";
        return 1;
    }

    results << "test_BFGS & True\n";
    return 0;
}

int test_BFGS2( std::ofstream &results ){
    /*!
     * Test of the BFGS optimization algorithm.
     *
     * :param std::ofstream &results: The output file
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>(lagrangian2);

    solverTools::floatVector x0 = { 0., 0., 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatOuts = { { .1, .2, .3, .4 } };
    intOuts = { { -1, -2 } };

    floatVector xAnswer = { -std::sqrt( 2. ) / 2, -std::sqrt( 2. ) / 2 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3}, {-0.4, -0.5, -0.6 }, { 7, 6, 5 } };
    intMatrix intOutsAnswer = { { -4 }, { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::BFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                        floatOuts, intOuts, floatArgs, intArgs,
                                        20, 1e-9, 1e-9, 1e-4, 5, true
                                      );

    if ( error ){
        error->print();
        results << "test_BFGS2 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( { x[ 0 ], x[ 1 ] }, xAnswer ) ){
        results << "test_BFGS2 (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) ){
        results << "test_BFGS2 (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) ){
        results << "test_BFGS2 (test 3) & False\n";
        return 1;
    }

    results << "test_BFGS2 & True\n";
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

    //Tests of the barrier homotopy solver
    test_aFxn( results );
    test_computeBarrierFunction( results );
    test_computeBarrierHomotopyResidual( results );
    test_computeBarrierHomotopyResidual2( results );
    test_barrierHomotopySolver( results );

    //Tests of the BFGS optimizer
    test_BFGS( results );
    test_BFGS2( results );

    //Close the results file
    results.close( );

    return 0;
}
