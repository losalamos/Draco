//----------------------------------*-C++-*----------------------------------//
// PreCond.cc
// Geoffrey M. Furnish
// Wed Nov 26 16:37:05 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/PreCond.hh"

#include "linalg/Banded_Matrix.hh"

using dsxx::Mat1;

//---------------------------------------------------------------------------//
// Solve QL * x = b for x.
//---------------------------------------------------------------------------//

template<class Solver>
void PreCond<Solver>::Left_PreCond( Mat1<T>& x, Mat1<T>&b )
{
    Banded_Matrix<double,7>& A = solver->get_A();

    switch(method) {

    case 0:
    // This is "null preconditioning".
        x = b;
        break;

    case 1:
    // Here is Jacobi preconditioning ("diagonal scaling").
        for( int i=0; i < ncp; i++ )
            x(i) = b(i) / A(i,3);
        break;

    default:
        throw "Unrecognized preconditioning option.";
    }
}

//---------------------------------------------------------------------------//
// Solve QR * x = b for x.
//---------------------------------------------------------------------------//

template<class Solver>
void PreCond<Solver>::Right_PreCond( Mat1<T>& x, Mat1<T>&b )
{
    x = b;
}

//---------------------------------------------------------------------------//
//                              end of PreCond.cc
//---------------------------------------------------------------------------//
