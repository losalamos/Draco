//----------------------------------*-C++-*----------------------------------//
// Input_edit.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Read, test, and manipulate input.
//---------------------------------------------------------------------------//

#include "sn/Input_edit.hh"

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

void Input_edit::read_edit_data()
{
    // instantiate an ifstream object and use it to read data

    ifstream input_file;

    input_file.open("snpp.in");

    if ( input_file.bad() )
    {
        cerr << "Error: Could not open snpp.in" << endl;
        exit (8);
    }

    input_file >> it_p   >> jt_p      >> kt_p    >> mm_p    >> isct_p;
    input_file >> dx_p   >> dy_p      >> dz_p    >> epsi_p;
    input_file >> ifxg_p >> iprint_p;
    input_file >> ibl_p  >> ibb_p     >> ibfr_p;

    input_file.close();

    // test for input errors and set parameters based on input

    if ( mm_p != 3 && mm_p != 6 )
    {
        cerr << "Error: mm must equal 3 or 6" << endl;
        exit (8);
    }

    if ( isct_p != 0 && isct_p != 1 )
    {
        cerr << "Error: isct must equal 0 or 1" << endl;
        exit (8);
    }

    isctp_p = isct_p + 1;
    nm_p    = isctp_p * isctp_p;
    maxop_p = 2 * mm_p * it_p;

    // Set boundary array dimensions based on boundary condition.

    if ( ibb_p == 1 )
        jbdim_p = 2 * mm_p * it_p;  // reflective boundary
    else
        jbdim_p = 1;                // vacuum boundary

    if ( ibfr_p == 1 )
        kbdim_p = 2 * mm_p * it_p;  // reflective boundary
    else
        kbdim_p = 1;                // vacuum boundary

    // Pre-calculate commonly used factors, for efficiency.

    itmm_p = it_p * mm_p;

    hi_p = new REAL [it_p];
    hj_p = new REAL [jt_p];
    hk_p = new REAL [kt_p];

    for ( int i=0 ; i < it_p ; i++ )
        hi_p[i] = 2.0 / dx_p;

    for ( int j=0 ; j < jt_p ; j++ )
        hj_p[j] = 2.0 / dy_p;

    for ( int k=0 ; k < kt_p ; k++ )
        hk_p[k] = 2.0 / dz_p;
}

Input_edit::~Input_edit()
{
    delete [] hi_p;
    delete [] hi_p;
    delete [] hi_p;
}

//---------------------------------------------------------------------------//
//                              end of Input_edit.cc
//---------------------------------------------------------------------------//

