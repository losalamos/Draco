//----------------------------------*-C++-*----------------------------------//
// snpp.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Main program for testing the structured sweeper.
//---------------------------------------------------------------------------//

// SN++ is a single node, ordered sweep, C++ version of:
//
// Method 3 Sweep - F90 version of SN++ within a PE. With blocks, resulting
//                  from spatial decompsition, being sweept by a diagonal
//                  plane, projected onto a 2D array of PE's. Serial in angle.

#include "sn/Input_edit.hh"
#include "sn/Inner_iter.hh"

#include <iostream.h>

int main()
{

    // set precision for output of floating point numbers

    cout.precision(13);

    // get input data from the input object

    Input_edit data;

    data.read_edit_data();

    // print output header to screen

    int isn;  // quadrature order, the N in SN

    if ( data.mm() == 3 )
        isn = 4;
    else
        isn = 6;

    cout << "Method 3 - Single Node C++ version" << endl;
    cout << " " << data.it() << " x " << data.jt() << " x " << data.kt() <<endl;
    cout << " S" << isn << "P" << data.isct() << endl;

    // do the inner iterations

    Inner_iter inner_iter_object;

    inner_iter_object.do_inner_iter( data );

    return (0);
}

//---------------------------------------------------------------------------//
//                              end of snpp.cc
//---------------------------------------------------------------------------//

