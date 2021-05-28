/**
  * \file test_solver_tools.cpp
  *
  * Tests for solver_tools
  */

#include<solver_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_solver_tools
#include <boost/test/included/unit_test.hpp>

typedef solverTools::errorOut errorOut;
typedef solverTools::errorNode errorNode;
typedef solverTools::floatType floatType;
typedef solverTools::floatVector floatVector;
typedef solverTools::floatMatrix floatMatrix;
typedef solverTools::intVector intVector;
typedef solverTools::intMatrix intMatrix;

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer )
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
     * A non-linear function for use in testing the solver. This function is a linear
     * function of the form
     * 
     * \f$ R = \left [ x + 1, y - 5.6 \right ]\f$
     * 
     * which has a solution at \f$\left ( -1, 5.6 \right )\f$
     * 
     * The function also sets floatOuts to
     * 
     * \f$\left [ \left [ -1 \right ], \left [ -1, -2, -3 \right ], \left [ 4, 5, 6 \right ] \right ]\f$
     * 
     * And intOuts to
     * 
     * \f$\left [ \left [ 1, 2, 8 \right ] \right ]\f$
     * 
     * \param &x: The variable vector. Of size 2.
     * \param &floatArgs: Floating point arguments to the function. None expected.
     * \param &intArgs: Integer arguments to the function. None expected.
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs.
     * \param &intOuts: Additional integer outputs.
     */

    if ( x.size( ) != 2 ){
        return new errorNode( "nlFnx1", "x must have a size of 2" );
    }

    floatType x0 = -1;
    floatType y0 = 5.6;

    residual = { x[ 0 ] - x0, x[ 1 ] - y0 };
    jacobian = { { 1, 0 }, { 0, 1 } };

    floatOuts = { { -1 }, { -1, -2, -3 }, { 4, 5, 6 } };
    intOuts = { { 1, 2, 8 } };

    return NULL;
}

errorOut nlFxn1( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                 floatVector &residual ){
    /*!
     * A non-linear function for use in testing the solver. An overload that hides the jacobian,
     * floatOut, and intOut arrays.
     * 
     * \param &x: The variable vector. Of size 2
     * \param &floatArgs: Floating point arguments to the function
     * \param &intArgs: Integer arguments to the function
     * \param &residual: The residual vector output.
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
     * A non-linear function for use in testing the solver. A polynomial function of the form
     * 
     * \f$ R = \left [ ( x - 1 ) ( x - 7 ) y, ( y - 1 ) ( x - 3 ) z, x y z \right ]\f$
     *
     * The function also sets floatOuts to
     * 
     * \f$\left [ \left [ -1 \right ], \left [ -1, -2, -3 \right ], \left [ 4, 5, 6 \right ] \right ]\f$
     * 
     * And intOuts to
     * 
     * \f$\left [ \left [ 1, 2, 8 \right ] \right ]\f$
     * 
     * \param &x: The variable vector. Of size 3.
     * \param &floatArgs: Floating point arguments to the function. Unused.
     * \param &intArgs: Integer arguments to the function. Unused.
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs.
     * \param &intOuts: Additional integer outputs.
     */

    if ( x.size( ) != 3 ){
        return new errorNode( "nlFxn2", "x must have a size of 3" );
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
     * A non-linear function for use in testing the solver. The same as the previously overloaded function
     * except this function obfuscates the computation of the Jacobian, the floatOuts, and the intOuts.
     * 
     * \param &x: The variable vector
     * \param &floatArgs: Floating point arguments to the function
     * \param &intArgs: Integer arguments to the function
     * \param &residual: The residual vector output.
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
     * require the use of the line-search algorithm. The function is of the form
     * 
     * \f$R = exp( -x ) - 1\f$
     * 
     * The function also sets floatOuts to
     * 
     * \f$\left [ \left [ -1 \right ], \left [ -1, -2, -3 \right ], \left [ 4, 5, 6 \right ] \right ]\f$
     * 
     * And intOuts to
     * 
     * \f$\left [ \left [ 1, 2 8 \right ] \right ]\f$
     * 
     * \param &x: The variable vector. Size 1.
     * \param &floatArgs: Floating point arguments to the function
     * \param &intArgs: Integer arguments to the function
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs.
     * \param &intOuts: Additional integer outputs.
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
     * require the use of the line-search algorithm. The same as the overloaded nlFxn3 except
     * this obfuscates the jacobian, the floatOuts, and the intOuts.
     * 
     * \param &x: The variable vector. Size 1.
     * \param &floatArgs: Floating point arguments to the function
     * \param &intArgs: Integer arguments to the function
     * \param &residual: The residual vector output.
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
     * require the use of the line-search algorithm. The function is of the form
     * 
     * \f$R = tanh( x )\f$
     *
     * \param &x: The variable vector. Of size 1.
     * \param &floatArgs: Floating point arguments to the function. Unused.
     * \param &intArgs: Integer arguments to the function. Unused.
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs. Unused.
     * \param &intOuts: Additional integer outputs. Unused.
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
     * the use of the bounded homotopy solver. This function performs error checking on the
     * floatArgs, intArgs, floatOuts, and intOuts arrays. The function is of the form.
     * 
     *  \f$ R = \left [ ( x - 1 ) ( x + 1 ) ( x + 1 ) ( x - 0.25 ) ( x + 0.1 ) \right ]\f$
     * 
     * The expected values for `floatArgs` are:
     * 
     * \f$ floatArgs = \left [ \left [ 0.1, 0.2, 0.3, 0.4 \right ], \left [ -0.01, -0.02 \right ] \right ] \f$
     * 
     * The expected values for `intArgs` are
     * 
     * \f$ intArgs = \left [ \left [ -1, -2, -3 \right ], \left [ 5, 4, 3, 2 \right ], \left [ 8, 9, 9 \right ] \right ]\f$
     * 
     * The expected incoming values for `floatOuts` are
     * 
     * \f$ floatOuts = \left [ \left [ 0, 1, 2 \right ], \left [ 7, -6 \right ], \left [ 0.24, 0.25 \right ] \right ]\f$
     *
     * The expected incoming values for `intOuts` are
     * 
     * \f$ intOuts = \left [ \left [ 1, 2, 3 \right ], \left [ -5, 6, 7, 8 \right ] \right ] \f$
     * 
     * The function currently throws an error if the expected values are not provided. `floatOuts` is 
     * updated to
     * 
     * \f$ floatOuts = \left [ \left [ 0.1, 1.1, 2.1 \right ], \left [ 7, -6 \right ], \left [ 0, 1, 2 \right ] \right ]\f$
     * 
     * `intOuts` is updated to
     * 
     * \f$ intOuts = \left [ \left [ -1, 0, 1 \right ], \left [ 1, 2, 3 \right ], \left [ -5, 6, 7, 8 \right ] \right ] \f$
     * 
     * \param &x: The variable vector. One value.
     * \param &floatArgs: Floating point arguments to the function
     * \param &intArgs: Integer arguments to the function
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs.
     * \param &intOuts: Additional integer outputs.
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
    if ( x.size( ) != 1 ){
        return new errorNode( "nlFxn5", "The x vector should have a size of 1" );
    }

    //floatArgs tests
    if ( floatArgs.size( ) != 2 ){
        return new errorNode( "nlFxn5", "The floatArgs matrix should have two values" );
    }

    if ( !vectorTools::fuzzyEquals( floatArgs[ 0 ], answer1 ) ){
        return new errorNode( "nlFxn5", "The first value of floatArgs should be { 0.1, 0.2, 0.3, 0.4 }" );
    }

    if ( !vectorTools::fuzzyEquals( floatArgs[ 1 ], answer2 ) ){
        return new errorNode( "nlFxn5", "The second value of floatArgs should be { -0.01, -0.02 }" );
    }

    //intArgs tests
    if ( intArgs.size( ) != 3 ){
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
    if ( floatOuts.size( ) != 3 ){
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
    if ( intOuts.size( ) != 2 ){
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
     * the use of the bounded homotopy solver. The function is of the form
     * 
     * \f$R = \left [ \left ( x - 1 \right ) \left ( x + 1 \right ) \left ( x - 0.25 \right ) \left ( x + 0.1 \right ), \left ( y - 1 \right ) \left ( y - 1 \right ), \left ( x + 5 \right ) \left ( z + 1 \right ) \right ]\f$
     *
     * \param &x: The variable vector. Three values required.
     * \param &floatArgs: Floating point arguments to the function. Unused.
     * \param &intArgs: Integer arguments to the function. Unused.
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs. Unused.
     * \param &intOuts: Additional integer outputs. Unused.
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
     * the use of the bounded homotopy solver. The function is of the form
     *
     * \f$R = log( x )\f$
     * 
     * \param &x: The variable vector. One value required.
     * \param &floatArgs: Floating point arguments to the function. Unused.
     * \param &intArgs: Integer arguments to the function. Unused.
     * \param &residual: The residual vector output.
     * \param &jacobian: The jacobian output.
     * \param &floatOuts: Additional floating point outputs. Unused.
     * \param &intOuts: Additional integer outputs. Unused.
     */

    if ( x.size( ) != 1 ){
        return new errorNode( "nlFxn7", "The x vector must have a size of 1" );
    }

    residual = { std::log( x[ 0 ] ) };

    jacobian = { { 1. / x[ 0 ] } };

    return NULL;
}

errorOut lagrangian1( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                      floatType &value, floatVector &gradient, floatMatrix &floatOuts, intMatrix &intOuts ){
    /*!
     * A lagrangian used to test the optimization tools. The function is of the form.
     * 
     * \f$ L = ( x - 1 ) ( x + 3 )\f$
     * 
     * `floatOuts` is updated to
     * 
     * \f$floatOuts = \left [ \left [ 1, 2, 3 \right ], \left [ -0.4, -0.5, -0.6 \right ] \right ] \f$
     * 
     * `intOuts` is updated to
     * 
     * \f$intOuts = \left [ \left [ 5, 6, 7 \right ], \left [ 8 \right ], \left [ 9, 10 \right ] \right ]\f$
     *
     * \param &x: A vector of the variable to be solved. One value required.
     * \param &floatArgs: Additional floating point arguments to residual. Unused.
     * \param &intArgs: Additional integer arguments to the residual. Unused.
     * \param &value: The value of the Lagrangian
     * \param &gradient: The gradient of the Lagrangian
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */

    if ( x.size( ) != 1 ){
        return new errorNode( "lagrangian1", "The x vector must have a size of 1" );
    }

    value = ( x[ 0 ] - 1 ) * ( x[ 0 ] + 3 );
    gradient = { ( x[ 0 ] + 3 ) + ( x[ 0 ] - 1 ) };

    floatOuts = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 } };
    intOuts = { { 5, 6, 7 }, { 8 }, { 9, 10 } };

    return NULL;
}

errorOut lagrangian2( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                      floatType &value, floatVector &gradient, floatMatrix &floatOuts, intMatrix &intOuts ){
    /*!
     * A lagrangian used to test the optimization tools
     *
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &value: The value of the Lagrangian
     * \param &gradient: The gradient of the Lagrangian
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */

    if ( x.size( ) != 3 ){
        return new errorNode( "lagrangian2", "The x vector must have a size of 3" );
    }

    if ( floatOuts.size( ) != 1 ){
        return new errorNode( "lagrangian2", "The floatOuts must have a size of 1" );
    }

    if ( intOuts.size( ) != 1 ){
        return new errorNode( "lagrangian2", "The intOuts must have a size of 1" );
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 0 ], { 0.1, 0.2, 0.3, 0.4 } ) ){
        return new errorNode( "lagrangian2", "The first value of the floatOuts is incorrect" );
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 0 ], { -1, -2 } ) ){
        return new errorNode( "lagrangian2", "The first value of the intOuts is incorrect" );
    }

    floatType _x = x[ 0 ];
    floatType _y = x[ 1 ];
    floatType _L = x[ 2 ];

    value = _x + _y + _L * ( _x * _x + _y * _y - 1 );

    gradient = { 1 + 2 * _L * _x,
                 1 + 2 * _L * _y,
                 _x * _x + _y * _y - 1 };

    floatOuts = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 }, { 7, 6, 5 } };
    intOuts = { { -4 }, { 5, 6, 7 }, { 8 }, { 9, 10 } };

    return NULL;
}

errorOut lagrangian3( const floatVector &x, const floatMatrix &floatArgs, const intMatrix &intArgs,
                      floatType &value, floatVector &gradient, floatMatrix &floatOuts, intMatrix &intOuts
                    ){
    /*!
     * A lagrangian used to test the optimization tools. The function is of the form
     * 
     * \f$L = x^2 y + z * \left ( x^2 y^2 - 3 \right ) \f$
     * 
     * `floatOuts` is expected to have an incoming value of
     * 
     * \f$ floatOuts = \left [ \left [ 0.1, 0.2, 0.3, 0.4 \right ] \right ] \f$
     * 
     * `intOuts` is expected to have an incoming value of
     * 
     * \f$ intOuts = \left [ \left [ 0 \right ], \left [ -1, -2 \right ] \right ] \f$
     * 
     * `floatOuts is updated to
     * 
     * \f$ floatOuts = \left [ \left [ 1, 2, 3 \right ], \left [ -0.4, -0.5, -0.6 \right ], \left [ 7, 6, 5 \right ] \right ] \f$
     * 
     * `intOuts` is updated to
     * 
     * \f$ intOuts = \left [ \left [ -4 \right ], \left [ 5, 6, 7 \right ], \left [ 8 \right ], \left [ 9, 10 \right ] \right ] \f$ 
     *
     * \param &x: A vector of the variable to be solved.
     * \param &floatArgs: Additional floating point arguments to residual
     * \param &intArgs: Additional integer arguments to the residual
     * \param &value: The value of the Lagrangian
     * \param &gradient: The gradient of the Lagrangian
     * \param &floatOuts: Additional floating point values to return.
     * \param &intOuts: Additional integer values to return.
     */

    if ( x.size( ) != 3 ){
        return new errorNode( "lagrangian3", "The x vector must have a size of 3" );
    }

    if ( floatOuts.size( ) != 1 ){
        return new errorNode( "lagrangian3", "The floatOuts must have a size of 1" );
    }

    if ( intOuts.size( ) != 1 ){
        return new errorNode( "lagrangian3", "The intOuts must have a size of 1" );
    }

    if ( !vectorTools::fuzzyEquals( floatOuts[ 0 ], { 0.1, 0.2, 0.3, 0.4 } ) ){
        return new errorNode( "lagrangian3", "The first value of the floatOuts is incorrect" );
    }

    if ( !vectorTools::fuzzyEquals( intOuts[ 0 ], { -1, -2 } ) ){
        return new errorNode( "lagrangian3", "The first value of the intOuts is incorrect" );
    }

    floatType _x = x[ 0 ];
    floatType _y = x[ 1 ];
    floatType _L = x[ 2 ];

    value = _x * _x * _y + _L * ( _x * _x + _y * _y - 3 );

    gradient = { 2 * _x * _y + 2 * _L * _x,
                 _x * _x + 2 * _L * _y,
                 _x * _x + _y * _y - 3 };

    floatOuts = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 }, { 7, 6, 5 } };
    intOuts = { { -4 }, { 5, 6, 7 }, { 8 }, { 9, 10 } };

    return NULL;
}

BOOST_AUTO_TEST_CASE( testCheckTolerance ){
    /*!
     * Test the tolerance checking function.
     */

    floatVector R   = {  1,   2, 3.00000, -4.0 };
    floatVector tol = { 1.5, 2.1, 3.00001,  4.1 };
    bool result;

    errorOut error = solverTools::checkTolerance( R, tol, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( result );

    tol[ 0 ] = .98;

    error = solverTools::checkTolerance( R, tol, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( ! result );

    tol[ 0 ] = 1.5;
    tol[ 3 ] = 3.8;
    error = solverTools::checkTolerance( R, tol, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( ! result );

    tol = { 1.6, 2.5 };

    error = solverTools::checkTolerance( R, tol, result );

    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testNewtonRaphson ){
    /*!
     * Tests of the Newton-Raphson solver
     * 
     * \param &results: The output file
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

    BOOST_CHECK( ! error );

    floatVector Rtmp;
    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    error = nlFxn1( x, { }, { }, Rtmp, Jtmp, fO, iO );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( Rtmp, { 0, 0 } ) );

    //The second test
    x0 = { 1, 1, 1 };
    floatOut.clear( );
    intOut.clear( );
    fO.clear( );
    iO.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn2 );
    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { } );

    BOOST_CHECK( ! error );

    error = nlFxn2( x, { }, { }, Rtmp, Jtmp, fO, iO );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( Rtmp, { 0, 0, 0 } ) );

    //The third test
    x0 = { 3 };
    floatOut.clear( );
    intOut.clear( );
    fO.clear( );
    iO.clear( );
    
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn3 );
    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { } );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, { 0 } ) );

    //The fourth test. Tests the bounded Newton method
    x0 = { 10. };
    floatOut.clear( );
    intOut.clear( );
    fO.clear( );
    iO.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn7 );

    solverTools::solverType linearSolver;
    floatMatrix J;

    intVector boundVariableIndices = { 0 };
    intVector boundSigns = { 0 };
    floatVector boundValues = { 1e-9 };
    floatMatrix Jexp = { { 1. } };

    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { }, linearSolver, J );

    BOOST_CHECK( ! error );
    BOOST_CHECK( vectorTools::fuzzyEquals( x, { 1. } ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( J, Jexp ) );

    //The fifth test: This test makes sure that when a Newton-Raphson iteration fails it
    //correctly returns a failure to converge.
    x0 = { -5. };
    floatOut.clear( );
    intOut.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn4 );

    error = solverTools::newtonRaphson( func, x0, x, converged, fatalError, floatOut, intOut, { }, { }, 5 );

    BOOST_CHECK( ! converged );

}

BOOST_AUTO_TEST_CASE( testFiniteDifference ){
    /*!
     * Test the finite difference jacobian calculator.
     */

    //The first test
    floatVector x0 = { 1.5, 6 };
    floatMatrix J;

    solverTools::stdFncNLF func;
    func = static_cast<solverTools::NonLinearFunction>( nlFxn1 );
    solverTools::finiteDifference( func, x0, J, { }, { } );

    floatVector Rtmp;
    floatMatrix result;
    floatMatrix floatOuts;
    intMatrix intOuts;
    nlFxn1( x0, { }, { }, Rtmp, result, floatOuts, intOuts );

    BOOST_CHECK( vectorTools::fuzzyEquals( J, result ) );

    //The second test
    x0 = { 1, 1, 1 };
    floatOuts.clear( );
    intOuts.clear( );
    func = static_cast<solverTools::NonLinearFunction>( nlFxn2 );
    solverTools::finiteDifference( func, x0, J, { }, { } );
    nlFxn2( x0, { }, { }, Rtmp, result, floatOuts, intOuts );

    BOOST_CHECK( vectorTools::fuzzyEquals( J, result ) );

}

BOOST_AUTO_TEST_CASE( testCheckJacobian ){
    /*!
     * Test the jacobian checking utility.
     */

    //The first test

    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>( nlFxn1 );
    bool isGood = false;
    floatType eps = 1e-6;
    floatType tolr = 1e-6;
    floatType tola = 1e-6;
    bool suppressOutput = true;

    floatVector x0 = { 0, 0 };
    errorOut error = solverTools::checkJacobian( func, x0, { }, { }, isGood, eps, tolr, tola, suppressOutput );

    BOOST_CHECK( ! error );

    BOOST_CHECK( isGood );

    //The second test
    solverTools::stdFncNLFJ badfunc;
    badfunc = [ & ]( const floatVector &x_, const floatMatrix &floatArgs_, const intMatrix &intArgs_,
                            floatVector &r, floatMatrix &j, floatMatrix &fO, intMatrix &iO ){
        errorOut e = func( x_, floatArgs_, intArgs_, r, j, fO, iO );
        j[ 0 ][ 1 ] = 0.1;
        return e;
    }; 

    error = solverTools::checkJacobian( badfunc, x0, { }, { }, isGood, eps, tolr, tola, suppressOutput );

    BOOST_CHECK( ! error );

    BOOST_CHECK( ! isGood );

}

BOOST_AUTO_TEST_CASE( testCheckLSCriteria ){
    /*!
     * Test the line search criteria
     */

    floatVector R  = { 1, 2, 3, 4, 5, 6 };
    floatVector Rp = { 2, 3, 4, 5, 6, 7 };
    bool result;

    solverTools::checkLSCriteria( R, Rp, result );

    BOOST_CHECK( result );

    R[ 0 ] = 100;

    solverTools::checkLSCriteria( R, Rp, result );

    BOOST_CHECK( ! result );

}

BOOST_AUTO_TEST_CASE( testHomotopySolver ){
    /*!
     * Test the Homotopy solver.
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

    BOOST_CHECK( ! error );

    floatVector Rtmp;
    floatMatrix Jtmp;
    floatMatrix fO;
    intMatrix iO;
    error = nlFxn1( x, { }, { }, Rtmp, Jtmp, fO, iO );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( Rtmp, { 0, 0 } ) );

    //The second test
    x0 = { 1, 1, 1 };
    floatOuts.clear( );
    intOuts.clear( );
    fO.clear( );
    iO.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn2 );
    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { } );

    BOOST_CHECK( ! error );

    error = nlFxn2( x, { }, { }, Rtmp, Jtmp, fO, iO );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( Rtmp, { 0, 0, 0 } ) );

    //The third test
    x0 = { 3 };
    floatOuts.clear( );
    intOuts.clear( );
    fO.clear( );
    iO.clear( );
    
    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn3 );
    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { },
                                         20, 1e-9, 1e-9, 1e-4, 5, 0.2, 0.01 );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, { 0 } ) );

    //The fourth test
    x0 = { 3 };
    floatOuts.clear( );
    intOuts.clear( );

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn4 );

    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { },
                                         20, 1e-9, 1e-9, 1e-4, 4, 1.0, 0.1 );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, { 0 } ) );

    //The fifth test ( hard bounds )
    x0 = { 10 };
    floatOuts.clear( );
    intOuts.clear( );

    solverTools::solverType linearSolver;
    floatMatrix J, Jexp;

    Jexp = { { 1 } };

    intVector variableIndices = { 0 };
    intVector barrierSigns = { 0 };
    floatVector barrierValues = { 1e-9 };

    func = static_cast< solverTools::NonLinearFunctionWithJacobian >( nlFxn7 );

    error = solverTools::homotopySolver( func, x0, x, converged, fatalErrorFlag, floatOuts, intOuts, { }, { },
                                         linearSolver, J, variableIndices, barrierSigns, barrierValues,
                                         false, 20, 1e-9, 1e-9, 1e-4, 4, 1.0, 0.1 );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, { 1 } ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( J, Jexp ) );

    
}

BOOST_AUTO_TEST_CASE( test_aFxn ){
    /*!
     * Test the computation of the "\f$a\f$" parameter in the Barrier function.
     */

    floatType pseudoT = .72;
    floatType logAfxn = 5.2;

    floatType answer  = 42.26671935907283;

    floatType result;

    errorOut error = solverTools::aFxn( pseudoT, logAfxn, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    floatType dadT;

    error = solverTools::aFxn( pseudoT, logAfxn, result, dadT );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    floatType eps = 1e-6;
    floatType delta = eps * pseudoT + eps;

    floatType aP, aM;

    error = solverTools::aFxn( pseudoT + delta, logAfxn, aP );

    BOOST_CHECK( ! error );

    error = solverTools::aFxn( pseudoT - delta, logAfxn, aM );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( aP - aM ) / ( 2 * delta ), dadT ) );

}

BOOST_AUTO_TEST_CASE( test_computeBarrierFunction ){
    /*!
     * Test the computation of the boundary function
     */

    floatType x        = 0.4;
    floatType pseudoT  = 0.25;
    floatType logAmax  = 5;
    floatType b        = 0.14;

    floatType negativeSignAnswer = -0.5964638357684787;
    floatType positiveSignAnswer = 1.4780926435784547;

    floatType result;

    errorOut error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, false, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, negativeSignAnswer ) );

    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, true, result );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, positiveSignAnswer ) );

    //Test the Jacobians
    floatType dbdx, dbdt;

    //Test the Jacobians when the sign is negative
    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, false, result, dbdx, dbdt );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, negativeSignAnswer ) );

    floatType eps = 1e-6;

    floatType dx = eps * fabs( x ) + eps;
    floatType dt = eps * fabs( pseudoT ) + eps;

    floatType bP, bM;

    error = solverTools::computeBarrierFunction( x + dx, pseudoT, logAmax, b, false, bP );

    BOOST_CHECK( ! error );

    error = solverTools::computeBarrierFunction( x - dx, pseudoT, logAmax, b, false, bM );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dx ), dbdx ) );

    error = solverTools::computeBarrierFunction( x, pseudoT + dt, logAmax, b, false, bP );

    BOOST_CHECK( ! error );

    error = solverTools::computeBarrierFunction( x, pseudoT - dt, logAmax, b, false, bM );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dt ), dbdt ) );

    //Test the Jacobians when the sign is positive
    error = solverTools::computeBarrierFunction( x, pseudoT, logAmax, b, true, result, dbdx, dbdt );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, positiveSignAnswer ) );

    error = solverTools::computeBarrierFunction( x + dx, pseudoT, logAmax, b, true, bP );

    BOOST_CHECK( ! error );

    error = solverTools::computeBarrierFunction( x - dx, pseudoT, logAmax, b, true, bM );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dx ), dbdx ) );

    error = solverTools::computeBarrierFunction( x, pseudoT + dt, logAmax, b, true, bP );

    BOOST_CHECK( ! error );

    error = solverTools::computeBarrierFunction( x, pseudoT - dt, logAmax, b, true, bM );

    BOOST_CHECK( ! error );

    BOOST_CHECK( vectorTools::fuzzyEquals( ( bP - bM ) / ( 2 * dt ), dbdt ) );

}

BOOST_AUTO_TEST_CASE( test_computeBarrierHomotopyResidual ){
    /*!
     * Test the computation of the barrier homotopy residual
     */


    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>( nlFxn5 );

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

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( residualAnswer, residualResult ) );

    //Check that the non-homotopy outputs are as expected.
    BOOST_CHECK( floatOuts.size( ) == 4 );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 1 ], floatOutsDefault[ 0 ] + 0.1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 2 ], floatOutsDefault[ 1 ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts[ 3 ], floatOutsDefault[ 0 ] ) );

    BOOST_CHECK( intOuts.size( ) == 3 );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts[ 0 ], intOutsDefault[ 0 ] - 2 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts[ 1 ], intOutsDefault[ 0 ] ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts[ 2 ], intOutsDefault[ 1 ] ) );

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

    BOOST_CHECK( ! error  );

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0 - dx, floatArgs, intArgs, rM, JM,
                                                         floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    floatVector gradCol = ( rP - rM ) / ( 2 * dx[ 0 ] );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, jacobian[ 0 ] ) );

    //test drdt
    eps = 1e-7;
    floatType dt = eps * floatArgsDefault[ 0 ][ 0 ] + eps;

    floatArgs = floatArgsDefault;
    floatArgs[ 0 ][ 0 ] += dt;

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0, floatArgs, intArgs, rP, JP,
                                                         floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    floatArgs = floatArgsDefault;
    floatArgs[ 0 ][ 0 ] -= dt;

    floatOuts = floatOutsDefault;
    intOuts   = intOutsDefault;

    error = solverTools::computeBarrierHomotopyResidual( func, x0, floatArgs, intArgs, rM, JM,
                                                         floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    gradCol = ( rP - rM ) / ( 2 * dt );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, floatOuts[ 0 ], 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( test_computeBarrierHomotopyResidual2 ){
    /*!
     * Test for the computation of the barrier homotopy residual.
     */


    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>( nlFxn6 );

    floatVector          x = { 0.15, 0.1, -1.2 };
    floatType   pseudoTime = 0.24;

    floatVector logAMaxVals = { 10, 6 };
    floatVector bvals       = { 0.1, -1.1 };

    intVector   variableIndices = { 0, 2 };
    intVector   residualIndices = { 2, 1 };
    intVector   signs           = { 0, 1 };

    floatMatrix floatArgs = { { pseudoTime }, bvals, logAMaxVals };
    intMatrix   intArgs   = { variableIndices, residualIndices, signs };

    floatMatrix floatOuts = { { } };
    intMatrix   intOuts   = { { } };

    floatVector residualAnswer = { 0.0244375 ,  0.53651154, -0.97499937 };

    floatVector residualResult;
    floatMatrix jacobian;

    errorOut error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, residualResult, jacobian,
                                                                  floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( residualResult, residualAnswer ) );

    //Tests of the Jacobians

    //Test the Jacobian w.r.t. x
    floatType eps = 1e-7;
    for ( unsigned int i = 0; i < x.size( ); i++ ){
        floatVector delta( x.size( ), 0 );
        delta[ i ] = eps * fabs( x[ i ] ) + eps;

        floatVector rP, rM;

        error = solverTools::computeBarrierHomotopyResidual( func, x + delta, floatArgs, intArgs, rP, jacobian,
                                                             floatOuts, intOuts );

        BOOST_CHECK( ! error  );

        error = solverTools::computeBarrierHomotopyResidual( func, x - delta, floatArgs, intArgs, rM, jacobian,
                                                             floatOuts, intOuts );

        BOOST_CHECK( ! error  );

        floatVector gradCol = ( rP - rM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size( ); j++ ){
            BOOST_CHECK( vectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) );
        }
    }

    //Test the Jacobian w.r.t. t
    floatType dt = eps * fabs( pseudoTime ) + eps;

    floatVector rP, rM;

    floatArgs[ 0 ][ 0 ] = pseudoTime + dt;

    error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, rP, jacobian,
                                                         floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    floatArgs[ 0 ][ 0 ] = pseudoTime - dt;

    error = solverTools::computeBarrierHomotopyResidual( func, x, floatArgs, intArgs, rM, jacobian,
                                                         floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    floatVector gradCol = ( rP - rM ) / ( 2 * dt );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradCol, floatOuts[ 0 ] ) );

}

BOOST_AUTO_TEST_CASE( test_barrierHomotopySolver ){
    /*!
     * Test the barrier Homotopy solver. This solver enables the addition of
     * bounds to a non-linear solve which can help prevent solutions from being
     * pulled into undesirable domains without having to resort to computing the
     * Hessian of the residual function as would be required for optimization
     * based techniques.
     */

    solverTools::stdFncNLFJ func;
    func = static_cast<solverTools::NonLinearFunctionWithJacobian>( nlFxn5 );

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

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    solverTools::solverType linearSolver;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                         barrierValues, logAMaxValues, floatArgs, intArgs,
                                                         implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                         floatOuts, intOuts, linearSolver, jacobian,
                                                         20, 1e-9, 1e-9, 1e-4, 5, true );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

    floatMatrix jacobianResult;

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = func( result, floatArgs, intArgs, residualResult, jacobianResult, floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( jacobian, jacobianResult ) );

    x0 = { 0. };
    implicitRefine = true;

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                barrierValues, logAMaxValues, floatArgs, intArgs,
                                                implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts,
                                                20, 1e-9, 1e-9, 1e-4, 5, true );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = solverTools::barrierHomotopySolver( func, dt, x0, variableIndices, residualIndices, barrierSigns,
                                                barrierValues, logAMaxValues, floatArgs, intArgs,
                                                implicitRefine, result, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts, linearSolver, jacobian,
                                                20, 1e-9, 1e-9, 1e-4, 5, true );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

    floatOuts = floatOutsDefault;
    intOuts = intOutsDefault;

    error = func( result, floatArgs, intArgs, residualResult, jacobianResult, floatOuts, intOuts );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( jacobian, jacobianResult ) );

}

BOOST_AUTO_TEST_CASE( test_applyBoundaryLimitation ){
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

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( x0 + dx, xAnswer1 ) );

    dx = dxDefault;
    x0[ 3 ] = -0.9;
    dx[ 3 ] = -0.3;

    error = solverTools::applyBoundaryLimitation( x0, variableIndices, barrierSigns, barrierValues, dx );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( x0 + dx, xAnswer2 ) );

    dx = dxDefault;
    error = solverTools::applyBoundaryLimitation( x0, variableIndices, barrierSigns, barrierValues, dx, 1e-9, 1e-9, true );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( x0 + dx, xAnswer3 ) );

}

BOOST_AUTO_TEST_CASE( test_BFGS ){
    /*!
     * Test of the BFGS optimization algorithm.
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>( lagrangian1 );

    solverTools::floatVector x0 = { 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatVector xAnswer = { -1 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 } };
    intMatrix intOutsAnswer = { { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::BFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                        floatOuts, intOuts, floatArgs, intArgs );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, xAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_BFGS2 ){
    /*!
     * Test of the BFGS optimization algorithm.
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>( lagrangian2 );

    solverTools::floatVector x0 = { 0., 0., 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatOuts = { { .1, .2, .3, .4 } };
    intOuts = { { -1, -2 } };

    floatVector xAnswer = { -std::sqrt( 2. ) / 2, -std::sqrt( 2. ) / 2 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 }, { 7, 6, 5 } };
    intMatrix intOutsAnswer = { { -4 }, { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::BFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                        floatOuts, intOuts, floatArgs, intArgs,
                                        20, 1e-9, 1e-9, 1e-4, 5, true
                                      );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( { x[ 0 ], x[ 1 ] }, xAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_homotopyBFGS ){
    /*!
     * Test of the homotopy BFGS optimization algorithm.
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>( lagrangian1 );

    solverTools::floatVector x0 = { 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatVector xAnswer = { -1 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 } };
    intMatrix intOutsAnswer = { { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::homotopyBFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts, floatArgs, intArgs,
                                                100, 1e-9, 1e-9, 1e-4, 10, 1.0, 0.1, true, 1.0
                                              );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( x, xAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_homotopyBFGS2 ){
    /*!
     * Test of the BFGS optimization algorithm.
     */

    solverTools::stdFncLagrangianG func;
    func = static_cast<solverTools::LagrangianFunctionWithGradient>( lagrangian2 );

    solverTools::floatVector x0 = { 0., 0., 0. };
    solverTools::floatVector x;

    bool convergeFlag, fatalErrorFlag;
    solverTools::floatMatrix floatArgs, floatOuts;
    solverTools::intMatrix intArgs, intOuts;

    floatOuts = { { .1, .2, .3, .4 } };
    intOuts = { { -1, -2 } };

    floatVector xAnswer = { -std::sqrt( 2. ) / 2, -std::sqrt( 2. ) / 2 };
    floatMatrix floatOutsAnswer = { { 1, 2, 3 }, { -0.4, -0.5, -0.6 }, { 7, 6, 5 } };
    intMatrix intOutsAnswer = { { -4 }, { 5, 6, 7 }, { 8 }, { 9, 10 } };

    errorOut error = solverTools::homotopyBFGS( func, x0, x, convergeFlag, fatalErrorFlag,
                                                floatOuts, intOuts, floatArgs, intArgs,
                                                100, 1e-9, 1e-9, 1e-4, 10, 1.0, 0.1, true, 1.0
                                              );

    BOOST_CHECK( ! error  );

    BOOST_CHECK( vectorTools::fuzzyEquals( { x[ 0 ], x[ 1 ] }, xAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( floatOuts, floatOutsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( intOuts, intOutsAnswer ) );

}
