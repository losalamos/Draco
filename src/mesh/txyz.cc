//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Geoffrey M. Furnish
// Wed May 13 09:53:44 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.cc"

#include "nml/Group.hh"

#include <iostream>
using namespace std;

int main( int argc, char *argv[] )
{
    cout << "t1: passed\n";

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    SP<Mesh_XYZ> spm = new Mesh_XYZ( mdb );

    Mesh_XYZ::cell_array x( spm ), y( spm ), z( spm );
    Mesh_XYZ::fcdsf xf( spm ), yf( spm ), zf( spm );

    x = 1.;
//     y = 2.;
    z = x + y;

    xf = 1.;
    yf = xf;
    zf = xf + yf;

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
