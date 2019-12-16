//Tests for constitutive_tools

#include<solver_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

typedef solverTools::errorOut errorOut;
typedef solverTools::floatType floatType;
typedef solverTools::floatVector floatVector;
typedef solverTools::floatMatrix floatMatrix;

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

    //Close the results file
    results.close();

    return 0;
}
